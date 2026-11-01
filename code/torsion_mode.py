"""Inter-ring torsional mode of a lanthanide double-decker from an ORCA Hessian (.hess).

usage: python torsion_mode.py name.hess [...]        e.g. calculations/freq/Lu_frq/Lu_frq.hess

Finds the normal mode with the largest overlap with a rigid counter-rotation of the two rings, and gives
its harmonic frequency, the reduced moment of inertia of the two rings, the root-mean-square torsional
amplitude at 298.15 K (quantum and classical harmonic oscillator) and the resulting change of the
thermally averaged angular factor (beta - cos 4 theta)/(1 + beta) at 45 deg (beta of the Ca form).
Reproduces the values quoted in Section S3 of the Supporting Information (LuPc2 23.7 cm-1, 2.7 deg, -1.5 %;
LaPc2 16.1 cm-1, 3.9 deg, -3.2 %)."""
import sys,numpy as np
AMU=1.66053907e-27; BOHR=0.529177; HB=1.054571817e-34; KB=1.380649e-23; C=2.99792458e10
def read_hess(f):
    L=open(f).read().split('\n'); i={l.strip():k for k,l in enumerate(L) if l.startswith('$')}
    k=i['$atoms']; n=int(L[k+1]); A=[L[k+2+j].split() for j in range(n)]
    el=[a[0] for a in A]; m=np.array([float(a[1]) for a in A]); X=np.array([[float(v) for v in a[2:5]] for a in A])*BOHR
    k=i['$vibrational_frequencies']; nf=int(L[k+1]); fr=np.array([float(L[k+2+j].split()[1]) for j in range(nf)])
    k=i['$normal_modes']; r,c=map(int,L[k+1].split()); M=np.zeros((r,c)); j=k+2
    while j<len(L) and not L[j].startswith('$') and L[j].strip():
        cols=list(map(int,L[j].split())); j+=1
        for q in range(r): v=L[j+q].split(); M[int(v[0]),cols]=[float(x) for x in v[1:]]
        j+=r
    return el,m,X,fr,M
def analyse(f,beta=0.1729,T=298.15):
    el,m,X,fr,M=read_hess(f); n=len(el); mt=[i for i,e in enumerate(el) if e not in ('H','C','N','O','S')][0]
    X=X-X[mt]; N8=[i for i in range(n) if el[i]=='N' and np.linalg.norm(X[i])<2.9]
    z=np.linalg.svd(X[[i for i in range(n) if el[i] in 'CN' and np.linalg.norm(X[i])<7]])[2][2]
    if np.mean(X[N8]@z)<0: z=-z
    s=np.sign(X@z); s[mt]=0
    t=np.array([s[i]*np.cross(z,X[i]) for i in range(n)])          # ring A turns +, ring B turns -
    w=(t*np.sqrt(m)[:,None]).ravel(); w/=np.linalg.norm(w)
    ov=[]
    for q in range(M.shape[1]):
        if fr[q]<1: continue
        v=(M[:,q].reshape(n,3)*np.sqrt(m)[:,None]).ravel(); nv=np.linalg.norm(v)
        if nv>0: ov.append((abs(w@v)/nv,fr[q],q))
    ov.sort(reverse=True)
    rp=np.linalg.norm(X-np.outer(X@z,z),axis=1)**2
    IA=sum(m[i]*rp[i] for i in range(n) if s[i]>0); IB=sum(m[i]*rp[i] for i in range(n) if s[i]<0); Ir=IA*IB/(IA+IB)*AMU*1e-20
    o,nu,q=ov[0]; om=2*np.pi*C*nu
    sq=HB/(2*Ir*om)/np.tanh(HB*om/(2*KB*T)); sc=KB*T/(Ir*om**2)
    for lab,s2 in (('quantum',sq),('classical',sc)):
        d=(np.exp(-8*s2)-1)/(1+beta)*100
        print('%s: mode %d, nu %.1f cm-1, overlap %.2f | I_red %.0f amu A2 | rms dtheta %.2f deg | <F>/F(45) - 1 = %+.2f %%'%(lab,q,nu,o,Ir/AMU/1e-20,np.degrees(np.sqrt(s2)),d))
    print('next candidates (nu, overlap):',' '.join('%.1f(%.2f)'%(b,a) for a,b,_ in ov[1:4]))
    print('imaginary frequencies:',int((fr<-1).sum()), 'lowest:',np.round(np.sort(fr[np.abs(fr)>1])[:4],1))
if __name__=='__main__':
    if len(sys.argv)<2: print(__doc__); sys.exit()
    for f in sys.argv[1:]: print('==',f); analyse(f)
