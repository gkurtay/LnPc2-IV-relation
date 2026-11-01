# Partition of NTO cube densities (hole, particle) onto metal, the two macrocycles and the substituents.
# usage: python nto_partition.py name_iv_hole.cube [...]   (expects name_iv_part.cube next to each hole cube)
# Output per structure: shares of the hole and particle densities and D_sub = substituent share (hole) - substituent share (particle).
# Grid points are assigned to the nearest atom (distance scaled by covalent radius); rings are the two largest
# fused ring systems of the bond graph (bridges removed); heteroatom = N, O, S outside the rings.
import sys,glob,os,numpy as np
R={1:.31,6:.76,7:.71,8:.66,9:.57,16:1.05,17:1.02,35:1.2}
def oku(f):
    L=open(f).read().split('\n'); t=L[2].split(); n=int(t[0]); o=np.array(t[1:4],float)
    ax=[L[3+i].split() for i in range(3)]; N=[int(a[0]) for a in ax]; V=np.array([a[1:4] for a in ax],float)
    A=np.array([L[6+i].split()[:5] for i in range(abs(n))],float); s=6+abs(n)+(1 if n<0 else 0)
    d=np.array(' '.join(L[s:]).split(),float)[:N[0]*N[1]*N[2]]
    return A[:,0].astype(int),A[:,2:5],o,N,V,d
def koprular(nb,dug):
    tin={};low={};B=set();c=[0]
    for s in dug:
        if s in tin: continue
        st=[(s,-1,iter(nb[s]))]; tin[s]=low[s]=c[0]; c[0]+=1
        while st:
            u,p,it=st[-1]; v=next(it,None)
            if v is None:
                st.pop()
                if st: w=st[-1][0]; low[w]=min(low[w],low[u]); B.add((w,u)) if low[u]>tin[w] else None
            elif v==p: continue
            elif v in tin: low[u]=min(low[u],tin[v])
            else: tin[v]=low[v]=c[0]; c[0]+=1; st.append((v,u,iter(nb[v])))
    return B
def parcala(Z,X):
    n=len(Z); M=[i for i in range(n) if Z[i]>20]; r=np.array([R.get(z,1.9) for z in Z])
    D=np.linalg.norm(X[:,None]-X[None],axis=2)
    nb={i:[j for j in range(n) if j!=i and j not in M and i not in M and D[i,j]<1.2*(r[i]+r[j])] for i in range(n)}
    g={i:set(nb[i]) for i in range(n) if i not in M}; ch=True
    while ch:
        ch=False
        for i in list(g):
            if len(g[i])<=1:
                for j in g[i]: g[j].discard(i)
                del g[i]; ch=True
    B=koprular({i:sorted(g[i]) for i in g},sorted(g)); B|={(b,a) for a,b in B}
    kal=set(g); bil=[]
    while kal:
        s=kal.pop(); k={s}; q=[s]
        while q:
            u=q.pop()
            for v in g[u]:
                if v in kal and (u,v) not in B: kal.discard(v); k.add(v); q.append(v)
        bil.append(k)
    bil.sort(key=len,reverse=True); c1,c2=bil[0],bil[1]
    ns=lambda c:sum(1 for i in c for j in nb[i] if j not in c1|c2 and Z[j]>1)
    if (ns(c1),X[list(c1)].sum(1).mean())<(ns(c2),X[list(c2)].sum(1).mean()): c1,c2=c2,c1
    et=np.full(n,3)
    for i in M: et[i]=0
    for i in c1: et[i]=1
    for i in c2: et[i]=2
    for i in range(n):
        if Z[i]==1 and et[i]==3:
            j=[k for k in nb[i]]; et[i]=et[j[0]] if j and et[j[0]] in (1,2) else 3
    return et,r
def kesir(f):
    Z,X,o,N,V,d=oku(f); et,r=parcala(Z,X*0.529177)
    g=np.stack(np.meshgrid(*[np.arange(k) for k in N],indexing='ij'),-1).reshape(-1,3)@V+o
    w=d**2; w/=w.sum(); top=np.zeros(4); het=0.
    for a in range(0,len(g),20000):
        p=g[a:a+20000]; k=np.argmin(np.linalg.norm(p[:,None]-X[None],axis=2)/r[None],axis=1)
        np.add.at(top,et[k],w[a:a+20000]); het+=w[a:a+20000][np.isin(Z[k],[7,8,16])&(et[k]==3)].sum()
    return top,het
def main(files):
    if not files:
        print('usage: python nto_partition.py name_iv_hole.cube [...]  (name_iv_part.cube must lie next to each hole cube)'); return
    for f in files:
        b=os.path.basename(f).replace('_iv_hole.cube',''); s=[]; sub={}
        for t in ('hole','part'):
            try:
                top,het=kesir(f.replace('_iv_hole','_iv_'+t)); sub[t]=top[3]
                s.append('%s M %.3f Pc1 %.3f Pc2 %.3f sub %.3f het %.3f'%(t,*top,het))
            except Exception as e: s.append('%s ERROR %s'%(t,e))
        ds=' | D_sub %+.3f'%(sub['hole']-sub['part']) if len(sub)==2 else ''
        print('%-24s %s%s'%(b,' | '.join(s),ds),flush=True)
if __name__=='__main__': main(sys.argv[1:])
