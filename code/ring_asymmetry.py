"""Ring spin asymmetry a = (s1 - s2)/(s1 + s2) of a lanthanide double-decker.

usage: python ring_asymmetry.py file [...]
  ORCA output (.out): Mulliken spin populations (last block), as used for the computed structures.
  Gaussian log (.log) with pop=hirshfeld: Hirshfeld spin densities, as used for the heteroleptic
  crystal structures of Section 3.9 (calculations/hirshfeld/).

s1 and s2 are the spin populations summed over the atoms on either side of the metal along the sandwich axis
(more than 0.3 A above or below it), i.e. each macrocycle with its substituents; the metal is excluded.
This is the definition behind ring_spin_asym_a in data/structures.csv and the values of the article.
The sign of a depends on which ring is taken as ring 1; the article reports |a|."""
import sys,os,re,numpy as np
H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,H)
import descriptors as D
SY=('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo '
    'Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu').split()
def gaussian(f):
    t=open(f,errors='ignore').read(); k='Standard orientation:' if 'Standard orientation:' in t else 'Input orientation:'
    E=[];X=[]
    for l in t[t.rfind(k):].split('\n')[5:]:
        p=l.split()
        if len(p)<6: break
        E.append(SY[int(p[1])-1]); X.append([float(v) for v in p[3:6]])
    k=t.rfind('Hirshfeld charges, spin densities')
    if k<0: raise SystemExit(f+': no Hirshfeld spin densities (run with pop=hirshfeld)')
    S=[float(l.split()[3]) for l in t[k:].split('\n')[2:2+len(E)]]
    return E,np.array(X),np.array(S),'Hirshfeld'
def orca(f):
    t=open(f,errors='ignore').read(); E,X=D.orca(f); k=t.rfind('MULLIKEN ATOMIC CHARGES AND SPIN')
    if k<0: raise SystemExit(f+': no Mulliken spin populations')
    S=[float(re.match(r'\s*\d+\s+\w+\s*:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',l).group(2)) for l in t[k:].split('\n')[2:2+len(E)]]
    return E,X,np.array(S),'Mulliken'
def asym(E,X,S):
    """published definition: spin populations summed on each side of the metal along the sandwich axis
    (atoms more than 0.3 A above or below the metal), a = (s1 - s2)/(s1 + s2)"""
    R=D.R; n=len(E); m=[i for i,e in enumerate(E) if e not in R][0]; X=X-X[m]
    N8=[i for i in range(n) if E[i]=='N' and np.linalg.norm(X[i])<2.9]
    z=np.linalg.svd(X[[i for i in range(n) if E[i] in 'CN' and np.linalg.norm(X[i])<7]])[2][2]
    top=[i for i in N8 if X[i]@z>0]; bot=[i for i in N8 if X[i]@z<0]
    ax=X[top].mean(0)-X[bot].mean(0); ax/=np.linalg.norm(ax); h=X@ax
    s=[S[h>0.3].sum(),S[h<-0.3].sum()]
    return (s[0]-s[1])/(s[0]+s[1]),s
if __name__=='__main__':
    if len(sys.argv)<2: print(__doc__); sys.exit()
    for f in sys.argv[1:]:
        E,X,S,kind=gaussian(f) if f.lower().endswith('.log') else orca(f)
        a,s=asym(E,X,S)
        print('%-28s %s  %s spin: ring 1 %.3f, ring 2 %.3f  |a| = %.3f'%(os.path.basename(f),E[[i for i,e in enumerate(E) if e not in D.R][0]],kind,s[0],s[1],abs(a)))
