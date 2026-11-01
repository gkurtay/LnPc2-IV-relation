"""Aryl substituents: torsion dependence of the IV band, its decomposition, and an additive per-aryl estimate.

usage:
  python aryl.py --series          reproduce the analysis of data/aryl_series.csv (Lu(Ph8Pc)2 at its X-ray geometry)
  python aryl.py structure.xyz     find aryl rings bonded to the benzo rings, report their torsions and the estimated shift

Analysis (series at fixed geometry, CAM-B3LYP/def2-TZVP): the shift of the IV band relative to the core
(all aryl -> H) is linear in the mean cos^2 of the aryl torsion phi. Writing the IV energy, a single beta HOMO->LUMO
excitation in all members, as the orbital gap plus an electron-hole term K (K = omega - gap) splits the shift into
an orbital part and a response part. The aryl share of the frontier orbitals (canonical cubes, same partition as
code/nto_partition.py) is compared with that of a methyl group at the same position, which sets the background.
Estimate for a new structure: shift(%) = sum_i [A0 + A2*cos^2(phi_i)] with A0, A2 per aryl from the 16-phenyl series."""
import sys,os,csv,math
import numpy as np
H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,H)
def fit(x,y):
    b,a=np.polyfit(x,y,1); return a,b,float(np.corrcoef(x,y)[0,1])
def series():
    R={r['id']:r for r in csv.DictReader(open(os.path.join(H,'..','data','aryl_series.csv')))}
    T=[R[k] for k in ('t35','xray','t70','t90')]; c0=R['core']; ev=lambda r:float(r['IV_cm-1'])/8065.544
    x=np.array([float(r['cos2_mean']) for r in T]); out={}
    for lab,f in (('IV shift (%)',lambda r:float(r['shift_vs_core_pct'])),('a1u splitting (%)',lambda r:100*(float(r['alpha_a1u_splitting_eV'])/float(c0['alpha_a1u_splitting_eV'])-1)),
                  ('IV shift (eV)',lambda r:ev(r)-ev(c0)),('gap change (eV)',lambda r:float(r['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV'])),
                  ('K change (eV)',lambda r:float(r['K_eV'])-float(c0['K_eV'])),('aryl share HOMO',lambda r:float(r['aryl_share_beta_HOMO']))):
        a,b,r=fit(x,[f(t) for t in T]); out[lab]=(a,b); print('%-18s = %+.4f %+.4f*cos2   r %.4f'%(lab,a,b,r))
    a,b=out['IV shift (eV)']; g=out['gap change (eV)'][0]; k=out['K change (eV)'][0]
    print('\ntorsion-independent part (cos2 = 0): %.4f eV = gap %.4f (%.0f%%) + K %.4f (%.0f%%)'%(a,g,100*g/a,k,100*k/a))
    m=R['me8']; dm=ev(m)-ev(c0); gm=float(m['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV'])
    print('methyl (same positions):            %.4f eV = gap %.4f (%.0f%%) + K %.4f'%(dm,gm,100*gm/dm,dm-gm))
    print('aryl share of beta HOMO/LUMO at 90 deg %s/%s, methyl %s/%s (partition background)'%(R['t90']['aryl_share_beta_HOMO'],R['t90']['aryl_share_beta_LUMO'],m['aryl_share_beta_HOMO'],m['aryl_share_beta_LUMO']))
    A0,A2=out['IV shift (%)'][0]/16,out['IV shift (%)'][1]/16
    p=R['ph4']; est=int(p['n_aryl'])*(A0+A2*float(p['cos2_mean']))
    print('\nper aryl: A0 %+.4f %%, A2 %+.4f %%;  Lu(Ph4Pc)2 estimate %+.2f %%, computed %s %%'%(A0,A2,est,p['shift_vs_core_pct']))
    return A0,A2
def torsions(E,X):
    r={'H':.31,'C':.76,'N':.71,'O':.66,'S':1.05}; n=len(E); d=np.linalg.norm(X[:,None]-X[None],axis=2)
    B=[[j for j in range(n) if j!=i and E[i] in r and E[j] in r and d[i,j]<1.2*(r[E[i]]+r[E[j]])] for i in range(n)]
    nH=lambda i:sum(E[j]=='H' for j in B[i]); nC=lambda i:sum(E[j]=='C' for j in B[i])
    cut={(i,j) for i in range(n) for j in B[i] if E[i]==E[j]=='C' and nH(i)==nH(j)==0 and nC(i)==nC(j)==3 and d[i,j]>1.45}
    S=[i for i in range(n) if E[i]=='N']; core=set(S)
    while S:
        i=S.pop()
        for j in B[i]:
            if (i,j) not in cut and j not in core: core.add(j); S.append(j)
    nv=lambda I:np.linalg.svd(X[I]-X[I].mean(0))[2][2]; c2=[]
    for i,j in [(i,j) for (i,j) in cut if i in core and j not in core]:
        ph={j}; st=[j]
        while st:
            a=st.pop()
            for b in B[a]:
                if b!=i and b not in ph and b not in core: ph.add(b); st.append(b)
        ring={i}
        for k in B[i]:
            if k in core and E[k]=='C': ring|={k}|{q for q in B[k] if q in core and E[q]=='C'}
        c2.append(float(np.dot(nv([k for k in ph if E[k]=='C']),nv(sorted(ring)))**2))
    return c2
if __name__=='__main__':
    a=sys.argv[1:]
    if not a or '--series' in a: series(); sys.exit()
    import io,contextlib; buf=io.StringIO()
    with contextlib.redirect_stdout(buf): A0,A2=series()
    import predict as P
    for f in a:
        E,X=P.D.read(f); o=P.D.desc(E,X); c2=torsions(E,X); pc,pn=P.rel(o['dN4'],o['dCa'],o['th'])
        s=sum(A0+A2*c for c in c2); tw=np.degrees(np.arccos(np.sqrt(np.clip(c2,0,1)))) if c2 else []
        print('%s  %d aryl rings  torsion %s  <cos2> %s'%(os.path.basename(f),len(c2),'%.1f-%.1f deg'%(min(tw),max(tw)) if c2 else '-','%.3f'%np.mean(c2) if c2 else '-'))
        print('  estimated aryl shift %+.2f %% | Ca form (core) %.0f cm-1 -> %.0f cm-1 (gas)'%(s,pc,pc*(1+s/100)))
