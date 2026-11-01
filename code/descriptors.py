# Geometric descriptors of a lanthanide double-decker (d_N4, d_Ca, theta), as used in the article.
# The same routine is applied to coordinates printed in ORCA outputs and to molecules cut from CIF files.
import numpy as np
from collections import deque
R={'H':.31,'C':.76,'N':.71,'O':.66,'S':1.05,'Cl':1.02,'F':.57,'Br':1.2}
def desc(E,X):
    n=len(E); M=[i for i,e in enumerate(E) if e not in R]; m=M[0]
    X=X-X[m]; r=np.linalg.norm(X,axis=1); d=np.linalg.norm(X[:,None]-X[None],axis=2)
    N8=[i for i in range(n) if E[i]=='N' and r[i]<2.9]
    z=np.linalg.svd(X[[i for i in range(n) if r[i]<7 and E[i] in 'CN']])[2][2]; pr=X[N8]@z
    top=[a for a,p in zip(N8,pr) if p>0]; bot=[a for a,p in zip(N8,pr) if p<0]
    assert len(top)==4 and len(bot)==4,(len(top),len(bot))
    def duz(N4):
        c=X[N4].mean(0); v=np.linalg.svd(X[N4]-c)[2][2]; return c,(v if v@c>0 else -v)
    (ct,nt),(cb,nb)=duz(top),duz(bot); ax=(ct-cb)/np.linalg.norm(ct-cb); dN4=(ct-cb)@ax
    e1=X[top[0]]-(X[top[0]]@ax)*ax; e1/=np.linalg.norm(e1); e2=np.cross(ax,e1); az=lambda i:np.degrees(np.arctan2(X[i]@e2,X[i]@e1))
    th1=np.mean([min(abs(((az(a)-az(q)+45)%90)-45) for q in bot) for a in top]); th2=np.mean([min(abs(((az(a)-az(q)+45)%90)-45) for a in top) for q in bot])
    th=(th1+th2)/2   # symmetric definition: mean of the two choices of reference ring
    h=[np.mean([(X[i]-c)@v for i in range(n) if E[i]=='C' and min(d[i,k] for k in N4)<1.45]) for N4,c,v in ((top,ct,nt),(bot,cb,nb))]
    return dict(M=E[m],dN4=float(dN4),dCa=float(dN4+h[0]+h[1]),th=float(th))
def molekul(f):
    import gemmi   # only needed for CIF input (pip install gemmi)
    doc=gemmi.cif.read(f); blk=[b for b in doc if len(b.find_loop('_atom_site_label'))][0]; st=gemmi.make_small_structure_from_block(blk); cell=st.cell
    sites=st.get_all_unit_cell_sites(); el=[s.element.name for s in sites]; fr=np.array([[s.fract.x,s.fract.y,s.fract.z] for s in sites])%1.0
    M=[i for i,e in enumerate(el) if e not in R]; m0=M[0]
    sh=np.array([[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)])
    F=(fr[None]+sh[:,None]).reshape(-1,3); E=el*len(sh)
    X=np.array([cell.orthogonalize(gemmi.Fractional(*p)).tolist() for p in F]); cm=np.array(cell.orthogonalize(gemmi.Fractional(*fr[m0])).tolist())
    s=np.linalg.norm(X-cm,axis=1)<16; X=X[s]; E=[E[i] for i in np.where(s)[0]]
    # remove duplicate sites
    keep=[]; 
    for i in range(len(E)):
        if all(np.linalg.norm(X[i]-X[j])>0.3 for j in keep): keep.append(i)
    X=X[keep]; E=[E[i] for i in keep]
    met=int(np.argmin(np.linalg.norm(X-cm,axis=1))); d=np.linalg.norm(X[:,None]-X[None],axis=2)
    N8=[i for i in range(len(E)) if E[i]=='N' and d[met,i]<2.9]
    B=[[j for j in range(len(E)) if j!=i and E[i] in R and E[j] in R and 0.6<d[i,j]<1.15*(R[E[i]]+R[E[j]])] for i in range(len(E))]
    bag=set(N8); q=deque(N8)
    while q:
        k=q.popleft()
        for j in B[k]:
            if j not in bag: bag.add(j); q.append(j)
    idx=[met]+sorted(bag); return [E[i] for i in idx],X[idx]
def xyz(f):
    L=open(f).read().split('\n'); n=int(L[0]); A=[l.split() for l in L[2:2+n]]
    return [a[0] for a in A],np.array([[float(v) for v in a[1:4]] for a in A])
def orca(f):
    """last CARTESIAN COORDINATES (ANGSTROEM) block of an ORCA output"""
    L=open(f,errors='ignore').read().split('\n'); B=[];cur=None
    for i,l in enumerate(L):
        if 'CARTESIAN COORDINATES (ANGSTROEM)' in l: cur=[]; B.append(cur); continue
        if cur is not None:
            t=l.split()
            if len(t)==4 and not t[0].startswith('-'): cur.append(t)
            elif cur and len(t)!=4: cur=None
    A=B[-1]; return [a[0] for a in A],np.array([[float(v) for v in a[1:4]] for a in A])
def read(f):
    if f.lower().endswith('.cif'): return molekul(f)
    if f.lower().endswith('.xyz'): return xyz(f)
    return orca(f)
