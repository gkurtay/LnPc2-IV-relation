"""Heteroleptic double-deckers: the relation with the ring-asymmetry and donor corrections, and the location of the hole.

usage:
  python heteroleptic.py structure.xyz|.out|.cif --a 0.20 [--dsub 0.005]   predict one structure
  python heteroleptic.py --table                                           reproduce the analysis of data/heteroleptic.csv

In a heteroleptic sandwich the two rings are inequivalent, the ground-state hole is shared unequally
(ring spin asymmetry a, e.g. from Hirshfeld spin populations) and, in the two-state picture, the band is raised
by 1/sqrt(1 - a^2). The prediction is  nu = nu_relation * (1 + c*D_sub) / sqrt(1 - a^2).
The ground-state hole sits on the ring with the larger share of the IV particle NTO (the IV transition moves
the electron into the hole); 'hole_on_donor_ring' in data/heteroleptic.csv is that share on the substituted ring.
Delta_eff = |a| * nu is the effective level offset between the two rings in the same picture.
Calibrated on computed geometries only; at the X-ray geometries of four measured heteroleptics the relation,
like TD-DFT itself, fails (see the article, Section 3.9)."""
import sys,os,csv,json,math
from collections import defaultdict
H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,H)
import predict as P
def corrected(d4,dc,t,a=0.0,dsub=0.0):
    pc,pn=P.rel(d4,dc,t,dsub); q=math.sqrt(1-a*a); return pc/q,pn/q
def table():
    R=list(csv.DictReader(open(os.path.join(H,'..','data','heteroleptic.csv'))))
    print('%-16s %-28s %6s %6s %6s %7s | dev Ca: raw  asym  both'%('id','donor class','a','Dsub','hole_D','Deff/eV'))
    for r in R:
        print('%-16s %-28s %+6.3f %6.3f %6.3f %7.3f |   %+5.1f %+5.1f %+5.1f'%(r['id'],r['donor_class'],float(r['ring_spin_asym_a']),float(r['D_sub']),
              float(r['hole_on_donor_ring']),float(r['Delta_eff_eV']),float(r['dev_Ca_pct']),float(r['dev_Ca_asym_corrected_pct']),float(r['dev_Ca_both_corrected_pct'])))
    G=defaultdict(list)
    for r in R: G[r['donor_class']].append(float(r['hole_on_donor_ring']))
    print('\nshare of the ground-state hole on the substituted ring, by donor class:')
    for k,v in sorted(G.items(),key=lambda kv:sum(kv[1])/len(kv[1])): print('  %-28s n %d  %.2f-%.2f'%(k,len(v),min(v),max(v)))
    for lab,k in (('raw','dev_Ca_pct'),('asymmetry corrected','dev_Ca_asym_corrected_pct'),('both corrections','dev_Ca_both_corrected_pct')):
        v=[float(r[k]) for r in R]; m=sum(v)/len(v); s=math.sqrt(sum((x-m)**2 for x in v)/len(v))
        print('Ca form, %-20s %+.1f to %+.1f %%  (mean %+.1f, sd %.1f)'%(lab,min(v),max(v),m,s))
if __name__=='__main__':
    a=sys.argv[1:]
    if '--table' in a or not a: table(); sys.exit()
    A=0.0; ds=0.0
    if '--a' in a: i=a.index('--a'); A=float(a[i+1]); del a[i:i+2]
    if '--dsub' in a: i=a.index('--dsub'); ds=float(a[i+1]); del a[i:i+2]
    for f in a:
        o=P.D.desc(*P.D.read(f)); pc,pn=corrected(o['dN4'],o['dCa'],o['th'],A,ds)
        print('%s  %s  d_N4 %.3f  d_Ca %.3f  theta %.1f  a %.3f  D_sub %.3f'%(os.path.basename(f),o['M'],o['dN4'],o['dCa'],o['th'],A,ds))
        for lab,p,k in (('Ca form',pc,'Ca_form'),('N4 form',pn,'N4_form')):
            s=P.C[k]['scale_exp']; print('  %-8s gas %6.0f cm-1 | solution (x%.3f) %6.0f cm-1'%(lab,p,s,p*s))
