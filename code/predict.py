"""Predict the IV band of a homoleptic lanthanide double-decker from its geometry.
usage: python predict.py structure.cif|structure.xyz|orca_output.out [--theta 45] [--dsub 0.0]
Prints d_N4, d_Ca, theta and the IV energy from both forms (gas phase, and scaled to solution).
--dsub: donor share of the IV transition (substituent share of the NTO hole minus that of the particle,
from code/nto_partition.py); 0 for structures without heteroatom donors (default)."""
import sys,os,json,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import descriptors as D
C=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data','constants.json')))
def rel(d4,dc,t,dsub=0.0):
    n,c=C['N4_form'],C['Ca_form']; a4=(n['beta']-np.cos(np.radians(4*t)))/(1+n['beta']); ac=(c['beta']-np.cos(np.radians(4*t)))/(1+c['beta'])
    return (c['A_cm-1']*np.exp(-c['kappa_per_A']*dc)*ac*(1+c['donor_c']*dsub), (n['p_cm-1']+n['q_cm-1_per_A']*d4)*a4*(1+n['donor_c']*dsub))
if __name__=='__main__':
    a=sys.argv[1:]; th=None; ds=0.0
    if '--theta' in a: i=a.index('--theta'); th=float(a[i+1]); del a[i:i+2]
    if '--dsub' in a: i=a.index('--dsub'); ds=float(a[i+1]); del a[i:i+2]
    for f in a:
        o=D.desc(*D.read(f)); t=o['th'] if th is None else th; pc,pn=rel(o['dN4'],o['dCa'],t,ds)
        print('%s  %s  d_N4 %.3f A  d_Ca %.3f A  theta %.1f deg%s  D_sub %.3f'%(os.path.basename(f),o['M'],o['dN4'],o['dCa'],o['th'],'' if th is None else ' (evaluated at %.1f)'%th,ds))
        for lab,p,k in (('Ca form',pc,'Ca_form'),('N4 form',pn,'N4_form')):
            s=C[k]['scale_exp']
            if t<C[k]['node_deg']: print('  %-8s below the coupling node (%.1f deg): the relation does not apply'%(lab,C[k]['node_deg'])); continue
            print('  %-8s gas %6.0f cm-1 | solution (x%.3f) %6.0f cm-1 = %5.0f nm'%(lab,p,s,p*s,1e7/(p*s)))
