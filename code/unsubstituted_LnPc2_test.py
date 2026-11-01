"""External test: the relation (constants and scale fixed) against the IV components of unsubstituted [Ln(Pc)2]
measured in CH2Cl2 by Ostendorp and Homborg, Z. Anorg. Allg. Chem. 1996, 622, 1222 (Table 2).
usage: python unsubstituted_LnPc2_test.py [--crystal-angles]"""
import sys,os,csv,json,numpy as np
H=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data')
C=json.load(open(os.path.join(H,'constants.json'))); n,c=C['N4_form'],C['Ca_form']
F=lambda b,t:(b-np.cos(np.radians(4*t)))/(1+b)
rN=lambda d,t:n['scale_exp']*(n['p_cm-1']+n['q_cm-1_per_A']*d)*F(n['beta'],t)
rC=lambda d,t:c['scale_exp']*c['A_cm-1']*np.exp(-c['kappa_per_A']*d)*F(c['beta'],t)
ang='--crystal-angles' in sys.argv; R={}
for r in csv.DictReader(open(os.path.join(H,'unsubstituted_LnPc2_test.csv'))):
    if 'not used' in r['structure']: continue
    t=float(r['theta_crystal_deg']) if ang else 45.0
    R.setdefault(r['metal'],dict(iv=[float(r[k]) for k in ('IV1_cm-1','IV2_cm-1')],N=[],C=[]))
    R[r['metal']]['N'].append(rN(float(r['d_N4_A']),t)); R[r['metal']]['C'].append(rC(float(r['d_Ca_A']),t))
dN=[];dC=[]
print('metal  n   N4 form  Ca form    IV1    IV2   dev N4 %  dev Ca %')
for M,v in R.items():
    pn,pc=np.mean(v['N']),np.mean(v['C']); e=v['iv'][0]; dN.append(100*(pn-e)/e); dC.append(100*(pc-e)/e)
    print('%-5s %2d  %7.0f  %7.0f  %5.0f  %5.0f   %+6.1f    %+6.1f'%(M,len(v['N']),pn,pc,e,v['iv'][1],dN[-1],dC[-1]))
print('mean absolute deviation from IV1: N4 %.1f %%, Ca %.1f %%  (theta = %s)'%(np.mean(np.abs(dN)),np.mean(np.abs(dC)),'crystal' if ang else '45 deg'))
