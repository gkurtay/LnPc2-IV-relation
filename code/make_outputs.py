"""Regenerate tables and figures of the article and its Supporting Information from the files in data/.

usage: python make_outputs.py [output_directory]          (default: ../outputs)
needs: numpy, scipy, matplotlib

Writes
  table1_constants.csv      Table 1: both forms fitted to the 11 donor-free structures, 95% bootstrap intervals, LOO errors
  fig2_rotation_scan.*      Figure 2: rigid rotation scan of LuPc2 with the Ca form and its bootstrap band
  fig3_loo.*                Figure 3: leave-one-out predictions of the fitting set
  fig6_aryl.*               Figure 6: aryl substituents (torsion dependence, gap/electron-hole decomposition, orbital shares)
  fig7_heteroleptic.*       Figure 7: heteroleptic double-deckers (location of the hole, asymmetry and donor corrections)
  tableS_unsubstituted_LnPc2.txt     Table S9c: external test against unsubstituted LnPc2 (Ostendorp and Homborg)
  tableS_aryl.txt, tableS_heteroleptic.txt   analyses printed by aryl.py and heteroleptic.py
  figS1-S3, figS5-S8        Supporting Information: bootstrap of the Ca form (S1), functional series (S2), solvent and
                            skew angle (S3), test against unsubstituted LnPc2 (S5), aryl series (S6), structure gallery (S7),
                            views along the molecular axis (S8). The data of Figure S4 (orbital scan) are in data/orbital_scan.csv.
Panel letters only; captions are in the article."""
import sys,os,csv,json,io,contextlib
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import font_manager
from scipy.optimize import least_squares
H=os.path.dirname(os.path.abspath(__file__)); DATA=os.path.join(H,'..','data'); sys.path.insert(0,H)
OUT=os.path.abspath(sys.argv[1] if len(sys.argv)>1 else os.path.join(H,'..','outputs')); os.makedirs(OUT,exist_ok=True)
INK='#1d1d1b'; INK2='#55554f'; MUT='#9b9a94'; HAIR='#d9d8d3'; WASH='#f1f0ec'
from figstyle import polish, MC as _MC, BG, STATBOX, chip
MC={k:_MC[k] for k in ('Lu','Y','La')}; MK={'Lu':'o','Y':'s','La':'^'}
FONT='Inter' if any('Inter' in f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family':FONT,'font.size':7,'axes.labelsize':7.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':0.7,
 'axes.edgecolor':INK2,'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'xtick.major.size':2.5,'ytick.major.size':2.5,
 'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'legend.fontsize':6.5,'savefig.dpi':600,'pdf.fonttype':42})
def harf(a,s,x=-0.13,y=1.0): a.text(x,y,s,transform=a.transAxes,fontsize=9,fontweight='bold',color=INK,va='bottom',ha='left')
def ygrid(a): a.grid(axis='y',color=HAIR,lw=0.45); a.set_axisbelow(True)
def save(fig,n):
    polish(fig)
    for e in ('png','pdf'): fig.savefig(os.path.join(OUT,n+'.'+e),bbox_inches='tight')
    plt.close(fig); print('wrote',n)
S={r['id']:r for r in csv.DictReader(open(os.path.join(DATA,'structures.csv')))}
FIT=['LuPc2bare_gas','YPc2bare_td','LaPc2bare_td','LuPc1_td','YPc1_td','LaPc1_td','YPc2_Me8_td','LuPc2bare_rot40','LuPc2bare_rot35','LuPc2bare_rot30','LuPc2bare_rot25']
assert set(FIT)=={k for k,r in S.items() if r['set']=='fitting set'}   # 11 donor-free structures; order fixes the bootstrap draws of the article
fl=lambda k,c:float(S[k][c])
t0=np.array([fl(k,'theta_deg') for k in FIT]); d4=np.array([fl(k,'d_N4_A') for k in FIT]); dc=np.array([fl(k,'d_Ca_A') for k in FIT]); nu=np.array([fl(k,'IV_cm-1') for k in FIT])
ang=lambda b,t:(b-np.cos(np.radians(4*np.asarray(t))))/(1+b)
gC=lambda q,d,t:q[0]*1e5*np.exp(-q[1]*d)*ang(q[2],t); gN=lambda q,d,t:(q[0]+q[1]*d)*ang(q[2],t)
FORMS={'Ca':(gC,dc,[4.67,1.29,0.17]),'N4':(gN,d4,[28000.,-6700.,0.16])}
def fitq(g,d,t,y,p0): return least_squares(lambda q:g(q,d,t)-y,p0).x
def table1(nboot=2000):
    rows=[]; P={}
    for f,(g,d,p0) in FORMS.items():
        rng=np.random.default_rng(7)
        p=fitq(g,d,t0,nu,p0); r=g(p,d,t0)-nu; P[f]=p
        loo=np.array([g(fitq(g,d[np.arange(len(nu))!=k],t0[np.arange(len(nu))!=k],nu[np.arange(len(nu))!=k],p),d[k],t0[k])-nu[k] for k in range(len(nu))])
        B=[]
        for _ in range(nboot):
            k=rng.choice(np.arange(len(nu)),len(nu))
            B.append(fitq(g,d[k],t0[k],nu[k],p))
        B=np.array(B); node=lambda b:np.degrees(np.arccos(np.clip(b,-1,1)))/4
        lo,hi=np.percentile(B,[2.5,97.5],axis=0); nlo,nhi=np.percentile(node(B[:,2]),[2.5,97.5])
        names=('A (1e5 cm-1)','kappa (1/A)','beta') if f=='Ca' else ('p (cm-1)','q (cm-1/A)','beta')
        for i,nm in enumerate(names): rows.append((f,nm,p[i],lo[i],hi[i]))
        rows+= [(f,'node (deg)',node(p[2]),nlo,nhi),(f,'RMS (cm-1)',np.sqrt(np.mean(r**2)),'',''),(f,'LOO RMS (cm-1)',np.sqrt(np.mean(loo**2)),'','')]
    with open(os.path.join(OUT,'table1_constants.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['form','quantity','value','ci_2.5','ci_97.5'])
        for x in rows: w.writerow([x[0],x[1]]+[('%.6g'%v if v!='' else '') for v in x[2:]])
    print('wrote table1_constants'); [print('  %-3s %-15s %12.5g  %s'%(x[0],x[1],x[2],'' if x[3]=='' else '[%.4g, %.4g]'%(x[3],x[4]))) for x in rows]
    return P
def fig2(P,nboot=1000):
    """Figure 2: rigid rotation scan of LuPc2: structures (Explorer), IV energy with the relation, oscillator strength, spin density at and above the node"""
    from matplotlib.gridspec import GridSpecFromSubplotSpec as SGS
    p=P['Ca']; node=np.degrees(np.arccos(p[2]))/4; DC=fl('LuPc2bare_gas','d_Ca_A'); rng=np.random.default_rng(7); tt=np.linspace(node,45,160); B=[]
    for _ in range(nboot):
        k=rng.choice(len(nu),len(nu)); q=fitq(gC,dc[k],t0[k],nu[k],p); B.append(gC(q,DC,tt))
    lo,hi=np.percentile(np.array(B),[2.5,97.5],axis=0)
    R=sorted([k for k in S if k.startswith('LuPc2bare_rot') and k!='LuPc2bare_rot45']+['LuPc2bare_gas'],key=lambda k:fl(k,'theta_deg'))
    th=np.array([fl(k,'theta_deg') for k in R]); v=np.array([fl(k,'IV_cm-1') for k in R]); f=np.array([fl(k,'f_osc') for k in R])
    a=np.abs(np.array([float(S[k]['ring_spin_asym_a'] or 'nan') for k in R]))
    fig=plt.figure(figsize=(7.0,5.0)); gs=GridSpec(3,3,figure=fig,width_ratios=[1,1,0.9],height_ratios=[0.95,1,1],wspace=0.38,hspace=0.42)
    # (a) structures along the axis
    g0=SGS(1,3,subplot_spec=gs[0,0:2],wspace=0.02)
    for k,(t,lab) in enumerate((('0','eclipsed'),('20','coupling node'),('45','staggered'))):
        ax=fig.add_subplot(g0[0,k]); off=img(ax,'rot_top_%s.png'%t); jf=os.path.join(EXI,'rot_top_%s.json'%t)
        if os.path.exists(jf) and off and t!='0': theta_arc(ax,'rot_top_%s.png'%t,json.load(open(jf)),off,lw=1.2)
        th_=json.load(open(jf))['geo']['th'] if os.path.exists(jf) else float(t)
        ax.text(0.5,-0.03,'θ = %.0f°, %s'%(abs(th_),lab),transform=ax.transAxes,ha='center',va='top',fontsize=6.6,color=INK)
        if k==0: harf(ax,'a',-0.05,0.98)
    # (b) IV energy
    A=fig.add_subplot(gs[1:,0:2])
    A.axvspan(24.5,46.5,color=MC['Lu'],alpha=0.05,lw=0,zorder=0); A.axvline(node,color=INK2,lw=0.6,ls=(0,(2.5,2)),zorder=1)
    A.text(35.5,11350,'fitted range',ha='center',fontsize=6.5,color=MC['Lu']); A.text(7.5,11350,'below the node',ha='center',fontsize=6.5,color=INK2); A.text(node+0.6,11350,'node %.1f°'%node,fontsize=6.5,color=INK2)
    A.fill_between(tt,lo,hi,color=MC['Lu'],alpha=0.16,lw=0,zorder=2); A.plot(tt,gC(p,DC,tt),color=MC['Lu'],lw=1.5,zorder=3)
    m=th>=24; A.scatter(th[m],v[m],s=28,c=MC['Lu'],edgecolors='white',linewidths=0.8,zorder=4)
    lo_=(th<24)&(np.abs(th-20)>0.5); A.scatter(th[lo_],v[lo_],s=26,facecolors='white',edgecolors=INK2,linewidths=0.9,zorder=4)
    k20=int(np.argmin(np.abs(th-20))); A.scatter([th[k20]],[v[k20]],s=30,marker='D',c=INK2,edgecolors='white',linewidths=0.8,zorder=4)
    A.annotate('ring-localised ground state\n|$a$| = %.2f, $f$ = %.3f'%(a[k20],f[k20]),xy=(th[k20],v[k20]),xytext=(1.0,700),fontsize=6.3,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.6,shrinkA=0,shrinkB=3))
    A.text(45.6,3300,'relation (Cα form)\nwith 95% bootstrap band',ha='right',fontsize=6.3,color=MC['Lu'])
    A.set_xlim(-1.5,46.5); A.set_xticks(range(0,46,5)); A.set_xlabel('skew angle θ (°)'); A.set_ylabel('IV energy (cm$^{-1}$)'); A.set_ylim(0,12000); ygrid(A); harf(A,'b',-0.1,1.02)
    # (c) oscillator strength
    Bf=fig.add_subplot(gs[0,2]); Bf.axvline(node,color=INK2,lw=0.6,ls=(0,(2.5,2)))
    Bf.plot(th,f,color=MC['Lu'],lw=1.1,zorder=2); Bf.scatter(th,f,s=12,c=MC['Lu'],edgecolors='white',linewidths=0.5,zorder=3)
    Bf.set_xlim(-1.5,46.5); Bf.set_xticks(range(0,46,15)); Bf.set_ylim(0,0.18); Bf.set_ylabel('oscillator strength $f$'); Bf.set_xlabel('θ (°)',labelpad=1); ygrid(Bf); harf(Bf,'c',-0.36,1.04)
    # (d) spin density at the node and above it
    g1=SGS(2,1,subplot_spec=gs[1:,2],hspace=0.25)
    for k,(t,kk) in enumerate((('20',k20),('25',int(np.argmin(np.abs(th-25)))))):
        ax=fig.add_subplot(g1[k,0]); img(ax,'rot_spin_%s.png'%t)
        ax.text(0.5,-0.02,'θ = %.1f°, |$a$| = %.2f'%(th[kk],a[kk]),transform=ax.transAxes,ha='center',va='top',fontsize=6.6,color=INK)
        if k==0: harf(ax,'d',-0.12,0.98)
    save(fig,'fig2_rotation_scan')
SHORT={'LuPc2bare_gas':'LuPc$_2$','YPc2bare_td':'YPc$_2$','LaPc2bare_td':'LaPc$_2$','LuPc1_td':'Lu-1','YPc1_td':'Y-1','LaPc1_td':'La-1','YPc2_Me8_td':'Y[Pc(Me)$_8$]$_2$',
       'LuPc2bare_rot40':'40°','LuPc2bare_rot35':'35°','LuPc2bare_rot30':'30°','LuPc2bare_rot25':'25°'}
def fig3(P):
    """Figure 3: leave-one-out test of both forms; inset, the cluster of optimised structures near 10000 cm-1"""
    from matplotlib.lines import Line2D
    fig=plt.figure(figsize=(7.0,3.7)); gs=GridSpec(2,2,figure=fig,height_ratios=[3.2,1],hspace=0.08,wspace=0.3); met=[S[k]['metal'] for k in FIT]
    rot=np.array(['rot' in k for k in FIT])
    for j,(fm,lab,s) in enumerate((('Ca','Cα form','a'),('N4','N$_4$ form','b'))):
        g,d,_=FORMS[fm]; L=np.array([g(fitq(g,d[np.arange(len(nu))!=k],t0[np.arange(len(nu))!=k],nu[np.arange(len(nu))!=k],P[fm]),d[k],t0[k]) for k in range(len(nu))])
        a=fig.add_subplot(gs[0,j]); r=fig.add_subplot(gs[1,j],sharex=a); lim=(2400,10900); e=np.sqrt(np.mean((L-nu)**2))
        def draw(ax,big=True):
            ax.fill_between(lim,[lim[0]-e,lim[1]-e],[lim[0]+e,lim[1]+e],color=WASH,lw=0,zorder=0); ax.plot(lim,lim,color=INK2,lw=0.6,ls=(0,(2.5,2)),zorder=1)
            for k in range(len(nu)):
                m=met[k]; kw=dict(s=34 if big else 30,marker=MK[m],zorder=3)
                if rot[k]: ax.scatter([nu[k]],[L[k]],facecolors='white',edgecolors=MC[m],linewidths=1.0,**kw)
                else: ax.scatter([nu[k]],[L[k]],c=MC[m],edgecolors='white',linewidths=0.7,**kw)
        draw(a)
        for k in range(len(nu)):
            if rot[k] and nu[k]<9000: a.annotate(SHORT[FIT[k]],(nu[k],L[k]),xytext=(6,-3),textcoords='offset points',fontsize=6,color=INK2,va='top')
            if FIT[k]=='LaPc2bare_td': a.annotate('LaPc$_2$, La-1',(nu[k],L[k]),xytext=(7,-3),textcoords='offset points',fontsize=6,color=INK2,ha='left',va='top')
        a.set_xlim(lim); a.set_ylim(lim); a.set_ylabel('LOO prediction (cm$^{-1}$)'); plt.setp(a.get_xticklabels(),visible=False); a.tick_params(axis='x',length=0)
        a.text(0.97,0.05,'%s\nLOO RMS %.0f cm$^{-1}$, n = %d'%(lab,e,len(nu)),transform=a.transAxes,ha='right',va='bottom',fontsize=6.6,color=INK,bbox=STATBOX,linespacing=1.4)
        # inset: the cluster near 10000 cm-1
        ia=a.inset_axes([0.03,0.54,0.42,0.44]); ia.set_facecolor('white'); ia.white=True; draw(ia,False); c0,c1=9230,10400; ia.set_xlim(c0,c1); ia.set_ylim(c0,c1)
        ia.set_xticks([]); ia.set_yticks([])
        for sp in ia.spines.values(): sp.set_visible(True); sp.set_color(MUT); sp.set_linewidth(0.5)
        ks=[k for k in range(len(nu)) if c0<nu[k]<c1]
        for side in (1,-1):   # 1: point below the diagonal, label to the right; -1: above, label to the left
            kk=sorted([k for k in ks if (1 if L[k]<=nu[k] else -1)==side],key=lambda k:L[k]); ys=[L[k] for k in kk]
            for it in range(50):
                for q in range(1,len(ys)):
                    if ys[q]-ys[q-1]<62: d_=(62-(ys[q]-ys[q-1]))/2; ys[q-1]-=d_; ys[q]+=d_
            for k,y in zip(kk,ys):
                ia.annotate(SHORT[FIT[k]],(nu[k],L[k]),xytext=(nu[k]+side*(70 if side>0 else 45),y),textcoords='data',fontsize=5.2,color=INK2,ha='left' if side>0 else 'right',va='center')
        a.indicate_inset_zoom(ia,edgecolor=MUT,alpha=0.6,lw=0.5)
        r.axhspan(-e,e,color=WASH,lw=0); r.axhline(0,color=INK2,lw=0.6,ls=(0,(2.5,2)))
        for k in range(len(nu)):
            m=met[k]; kw=dict(s=26,marker=MK[m],zorder=3)
            if rot[k]: r.scatter([nu[k]],[L[k]-nu[k]],facecolors='white',edgecolors=MC[m],linewidths=1.0,**kw)
            else: r.scatter([nu[k]],[L[k]-nu[k]],c=MC[m],edgecolors='white',linewidths=0.7,**kw)
        r.set_ylim(-520,520); r.set_yticks([-400,0,400]); r.set_ylabel('Δ (cm$^{-1}$)'); r.set_xlabel('computed IV energy (cm$^{-1}$)'); a.spines['bottom'].set_visible(False); harf(a,s,-0.17,1.02)
    H=[Line2D([],[],marker=MK[m],ls='',mfc=MC[m],mec='white',ms=5.5,label=m+'Pc$_2$ family, optimised') for m in ('Lu','Y','La')]+[Line2D([],[],marker='o',ls='',mfc='white',mec=MC['Lu'],mew=1.0,ms=5,label='LuPc$_2$, rigid rotation')]
    fig.legend(handles=H,loc='upper center',ncol=4,bbox_to_anchor=(0.52,1.03),fontsize=6.3,columnspacing=1.2,handletextpad=0.3)
    save(fig,'fig3_loo')
# ---------------- Supporting Information figures ----------------
def figS1(P,nboot=2000):
    rng=np.random.default_rng(7); B=np.array([fitq(gC,dc[k],t0[k],nu[k],P['Ca']) for k in (rng.choice(np.arange(len(nu)),len(nu)) for _ in range(nboot))])
    nd=np.degrees(np.arccos(np.clip(B[:,2],-1,1)))/4; q=lambda v:np.percentile(v,[0.5,99.5]); pad=lambda a,b:(a-0.08*(b-a),b+0.08*(b-a))
    fig=plt.figure(figsize=(7.0,2.4)); gs=GridSpec(1,3,figure=fig,wspace=0.45); A,Bx,Cx=[fig.add_subplot(gs[0,i]) for i in range(3)]
    A.scatter(B[:,1],B[:,0],s=3,c=MC['Lu'],alpha=0.25,lw=0); A.scatter([P['Ca'][1]],[P['Ca'][0]],s=30,c=INK,edgecolors='white',zorder=3)
    A.set_xlim(*pad(*q(B[:,1]))); A.set_ylim(*pad(*q(B[:,0]))); A.set_xlabel('κ (Å$^{-1}$)'); A.set_ylabel('A (10$^5$ cm$^{-1}$)'); A.text(0.97,0.05,('r = %.3f'%np.corrcoef(B[:,0],B[:,1])[0,1]).replace('-','−'),transform=A.transAxes,ha='right',fontsize=6.5,color=INK2,bbox=STATBOX); harf(A,'a',-0.32,1.02)
    Bx.scatter(B[:,1],B[:,2],s=3,c=MC['Lu'],alpha=0.25,lw=0); Bx.scatter([P['Ca'][1]],[P['Ca'][2]],s=30,c=INK,edgecolors='white',zorder=3)
    Bx.set_xlim(*pad(*q(B[:,1]))); Bx.set_ylim(*pad(*q(B[:,2]))); Bx.set_xlabel('κ (Å$^{-1}$)'); Bx.set_ylabel('β'); Bx.text(0.97,0.05,('r = %.3f'%np.corrcoef(B[:,1],B[:,2])[0,1]).replace('-','−'),transform=Bx.transAxes,ha='right',fontsize=6.5,color=INK2,bbox=STATBOX); harf(Bx,'b',-0.32,1.02)
    a_,b_=q(nd); Cx.hist(nd[(nd>=a_)&(nd<=b_)],bins=30,color=MC['Lu'],alpha=0.8); lo,hi=np.percentile(nd,[2.5,97.5]); [Cx.axvline(v,color=INK2,lw=0.7,ls=(0,(2.5,2))) for v in (lo,hi)]
    Cx.set_xlabel('coupling node θ$_0$ (°)'); Cx.set_ylabel('bootstrap resamples'); ygrid(Cx); harf(Cx,'c',-0.32,1.02)
    save(fig,'figS1_bootstrap')
FUN=[('CAM-B3LYP',INK),('B3LYP',MC['Lu']),('PBE0',MC['Y']),('ωB97X-D3',MC['La']),('BLYP35','#b0417a')]
def figS2():
    ids={'CAM-B3LYP':['LuPc2bare_rot25','LuPc2bare_rot30','LuPc2bare_rot35','LuPc2bare_rot40','LuPc2bare_gas']}
    for f,suf in (('B3LYP','b3lyp'),('PBE0','pbe0'),('ωB97X-D3','wb97xd'),('BLYP35','blyp35')): ids[f]=['LuRot%d_%s'%(t,suf) for t in (25,30,35,40)]+['LuPc2bare_'+suf]
    fig=plt.figure(figsize=(7.0,2.6)); gs=GridSpec(1,2,figure=fig,wspace=0.32); A,Bx=fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]); tt=np.linspace(15,45,120); rows=[]
    for f,c in FUN:
        th=np.array([fl(k,'theta_deg') for k in ids[f]]); v=np.array([fl(k,'IV_cm-1') for k in ids[f]]); y=v/v[-1]
        b=least_squares(lambda q:ang(q[0],th)-y,[0.2]).x[0]; nd=np.degrees(np.arccos(b))/4; rows.append((f,b,nd))
        A.plot(tt,ang(b,tt),color=c,lw=1.1); A.scatter(th,y,s=18,c=c,edgecolors='white',linewidths=0.6,zorder=3,label='%s (node %.1f°)'%(f,nd))
        Bx.scatter(th,v,s=18,c=c,edgecolors='white',linewidths=0.6,zorder=3); Bx.plot(th,v,color=c,lw=0.8,alpha=0.6)
    A.axhline(0,color=HAIR,lw=0.7); A.set_xlim(15,46); A.set_xlabel('skew angle θ (°)'); A.set_ylabel('ν(θ) / ν(45°)'); A.legend(loc='lower right',fontsize=6.0); ygrid(A); harf(A,'a',-0.2,1.02)
    Bx.set_xlim(22,46); Bx.set_xlabel('skew angle θ (°)'); Bx.set_ylabel('IV energy (cm$^{-1}$)'); ygrid(Bx); harf(Bx,'b',-0.2,1.02)
    save(fig,'figS2_functionals'); return rows
def figS3():
    R=[r for r in csv.DictReader(open(os.path.join(DATA,'solvent.csv'))) if r.get('solvent','chlorobenzene')=='chlorobenzene']   # Figure S3: chlorobenzene series
    fig,A=plt.subplots(figsize=(3.4,2.5))
    for m in ('Lu','La'):
        rr=[r for r in R if r['metal']==m]; A.scatter([float(r['theta_deg']) for r in rr],[float(r['ratio']) for r in rr],s=28,marker=MK[m],c=MC[m],edgecolors='white',linewidths=0.7,zorder=3,label=m)
    b=[r for r in R if r['id_solvent'] in ('LuPc2bare_cpcm','LuPc2bare_r30_cb_td2')]; b.sort(key=lambda r:float(r['theta_deg']))
    A.annotate('',xy=(float(b[0]['theta_deg']),float(b[0]['ratio'])),xytext=(float(b[1]['theta_deg']),float(b[1]['ratio'])),arrowprops=dict(arrowstyle='->',color=INK2,lw=0.8,shrinkA=5,shrinkB=5))
    A.text(38.5,0.903,'LuPc$_2$ rotated 45° → 30°',fontsize=6.0,color=INK2,ha='center')
    A.set_xlim(27,47); A.set_ylim(0.86,0.96); A.set_xlabel('skew angle θ (°)'); A.set_ylabel('ν(CPCM, chlorobenzene) / ν(gas)'); A.legend(loc='lower right'); ygrid(A)
    save(fig,'figS3_solvent')
LNO=['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
def figS4():
    C=json.load(open(os.path.join(DATA,'constants.json'))); c,n=C['Ca_form'],C['N4_form']; M={}
    for r in csv.DictReader(open(os.path.join(DATA,'unsubstituted_LnPc2_test.csv'))):
        m=M.setdefault(r['metal'],dict(iv=(float(r['IV1_cm-1']),float(r['IV2_cm-1'])),C=[],N=[]))
        if 'not used' in r['structure']: continue
        m['C'].append(c['scale_exp']*c['A_cm-1']*np.exp(-c['kappa_per_A']*float(r['d_Ca_A']))); m['N'].append(n['scale_exp']*(n['p_cm-1']+n['q_cm-1_per_A']*float(r['d_N4_A'])))
    L=[m for m in LNO if m in M]; x=np.arange(len(L))
    fig,A=plt.subplots(figsize=(4.6,2.7))
    A.plot(x,[M[m]['iv'][0] for m in L],color=INK,lw=1.0); A.scatter(x,[M[m]['iv'][0] for m in L],s=24,c=INK,edgecolors='white',zorder=3,label='IV1, lowest component')
    A.plot(x,[M[m]['iv'][1] for m in L],color=MUT,lw=1.0); A.scatter(x,[M[m]['iv'][1] for m in L],s=24,facecolors='white',edgecolors=MUT,linewidths=1.0,zorder=3,label='IV2')
    for k,col,mk,lab in (('C',MC['Lu'],'D','Cα form × scale (θ = 45°)'),('N',MC['Y'],'s','N$_4$ form × scale (θ = 45°)')):
        xx=[i for i,m in zip(x,L) if M[m][k]]; A.scatter(xx,[np.mean(M[m][k]) for m in L if M[m][k]],s=22,marker=mk,c=col,edgecolors='white',linewidths=0.6,zorder=4,label=lab)
    A.set_xticks(x); A.set_xticklabels(L); A.set_ylabel('IV energy (cm$^{-1}$)'); A.legend(loc='lower right',fontsize=6.0); ygrid(A)
    save(fig,'figS5_unsubstituted_LnPc2')
def figS5():
    R={r['id']:r for r in csv.DictReader(open(os.path.join(DATA,'aryl_series.csv')))}; T=[R[k] for k in ('t35','xray','t70','t90')]; c0=R['core']
    x=np.array([float(r['cos2_mean']) for r in T]); xx=np.linspace(0,0.75,50)
    fig=plt.figure(figsize=(7.0,2.5)); gs=GridSpec(1,3,figure=fig,wspace=0.5); A,Bx,Cx=[fig.add_subplot(gs[0,i]) for i in range(3)]
    for ax,k,lab,ref in ((A,'alpha_a1u_splitting_eV','α a$_{1u}$ splitting, change (%)',float(c0['alpha_a1u_splitting_eV'])),(Bx,'beta_HOMO_LUMO_gap_eV','β HOMO–LUMO gap, change (%)',float(c0['beta_HOMO_LUMO_gap_eV']))):
        y=np.array([100*(float(r[k])/ref-1) for r in T]); b,a=np.polyfit(x,y,1); ax.plot(xx,a+b*xx,color=INK,lw=1.0); ax.scatter(x,y,s=22,c=INK,edgecolors='white',zorder=3,label='16 phenyl, rotated')
        for rid,mk,lb,xv in (('me8','s','methyl',0.0),('ph4','D','8 phenyl (X-ray torsions)',float(R['ph4']['cos2_mean']))):
            ax.scatter([xv],[100*(float(R[rid][k])/ref-1)],s=24,marker=mk,facecolors='white',edgecolors=INK2,linewidths=1.0,zorder=3,label=lb)
        ax.axhline(0,color=HAIR,lw=0.7); ax.set_xlabel('mean cos²φ'); ax.set_ylabel(lab); ygrid(ax)
    A.legend(loc='lower left',fontsize=5.8); harf(A,'a',-0.32,1.02); harf(Bx,'b',-0.32,1.02)
    f=np.array([float(r['f_osc']) for r in T]); b,a=np.polyfit(x,f,1); Cx.plot(xx,a+b*xx,color=INK,lw=1.0); Cx.scatter(x,f,s=22,c=INK,edgecolors='white',zorder=3)
    Cx.scatter([0],[float(R['me8']['f_osc'])],s=24,marker='s',facecolors='white',edgecolors=INK2,linewidths=1.0,zorder=3); Cx.axhline(float(c0['f_osc']),color=INK2,lw=0.7,ls=(0,(2.5,2)))
    Cx.text(0.74,float(c0['f_osc'])+0.0015,'core',ha='right',fontsize=6.0,color=INK2); Cx.set_xlabel('mean cos²φ'); Cx.set_ylabel('oscillator strength $f$ of the IV band'); ygrid(Cx); harf(Cx,'c',-0.32,1.02)
    save(fig,'figS6_aryl_orbitals')

EXI=os.path.join(OUT,'explorer')   # structure panels exported from the LnPc2 Explorer (code/explorer_images.js)
def img(ax,name,box=None):
    """draw an Explorer export, cropped to its content; box = (x0,x1,y0,y1) fractions of the cropped image"""
    import matplotlib.image as mpimg
    f=os.path.join(EXI,name)
    ax.axis('off')
    if not os.path.exists(f): ax.text(0.5,0.5,'(run code/explorer_images.js)',ha='center',va='center',fontsize=6,color=MUT,transform=ax.transAxes); return
    a=mpimg.imread(f)[...,:3]; m=a.mean(2)<0.985; ys,xs=np.where(m); a=a[ys.min():ys.max()+1,xs.min():xs.max()+1]
    if box: h,w=a.shape[:2]; a=a[int(box[2]*h):int(box[3]*h),int(box[0]*w):int(box[1]*w)]
    ax.imshow(a,interpolation='lanczos'); ax.set_aspect('equal')
    return xs.min()+(int(box[0]*(xs.max()-xs.min()+1)) if box else 0), ys.min()+(int(box[2]*(ys.max()-ys.min()+1)) if box else 0)
def ilab(ax,t,y=1.0): ax.text(0.5,y,t,transform=ax.transAxes,ha='center',va='bottom' if y>0.5 else 'top',fontsize=7,color=INK)

def metal_fix(name,scr,off):
    """3Dmol screen projection carries a container offset; anchor it on the metal sphere (Explorer colour #2a9d6f) found in the image"""
    import matplotlib.image as mpimg
    a=mpimg.imread(os.path.join(EXI,name))[...,:3]; t=np.array([0x2a,0x9d,0x6f])/255.
    m=np.linalg.norm(a-t,axis=2)<0.12; ys,xs=np.where(m)
    if len(xs)<20: return np.array(scr)-np.array(off)
    c=np.array([np.median(xs),np.median(ys)]); S=np.array(scr,float); return S-S[0]+c-np.array(off)

def theta_arc(ax,name,js,off,lw=1.6,fs=7):
    S=metal_fix(name,js['scr'],off); m=S[0]; t=S[1]; b=min(S[5:9],key=lambda q:np.linalg.norm((q-m)/np.linalg.norm(q-m)-(t-m)/np.linalg.norm(t-m)))
    for q,c in ((t,'#1d3f9e'),(b,'#a8541a')): ax.plot([m[0],m[0]+3.7*(q[0]-m[0])],[m[1],m[1]+3.7*(q[1]-m[1])],color=c,lw=lw,solid_capstyle='round')
    a0=np.degrees(np.arctan2(-(t[1]-m[1]),t[0]-m[0])); a1=np.degrees(np.arctan2(-(b[1]-m[1]),b[0]-m[0])); R=3.4*np.linalg.norm(t-m); lo,hi=sorted([a0,a1])
    if hi-lo>180: lo,hi=hi,lo+360
    tt=np.radians(np.linspace(lo,hi,40)); ax.plot(m[0]+R*np.cos(tt),m[1]-R*np.sin(tt),color=INK,lw=1.0)
def figS6():
    L={k:r['label'] for k,r in S.items()}
    G=['LuPc2bare_gas','YPc2bare_td','LaPc2bare_td','LuPc1_td','LuPc2_td','LuPc3_td2','LuPc4_td2','LuPc8_td','LaPc2_td']
    fig=plt.figure(figsize=(7.0,7.2)); gs=GridSpec(3,3,figure=fig,hspace=0.32,wspace=0.06)
    for k,i in enumerate(G):
        ax=fig.add_subplot(gs[k//3,k%3]); img(ax,'gal_S6_%s.png'%i); f=os.path.join(EXI,'gal_S6_%s.json'%i)
        if os.path.exists(f):
            j=json.load(open(f)); g=j['geo']; lnn=[x[2] for x in j['m']]
            ax.text(0.5,-0.02,'%s\nθ %.1f°, $d_{\\mathrm{N4}}$ %.3f Å, $d_{\\mathrm{C\\alpha}}$ %.3f Å'%(L.get(i,i),g['th'],g['dN4'],g['dCa']),transform=ax.transAxes,ha='center',va='top',fontsize=6.4,color=INK)
        harf(ax,'abcdefghi'[k],0.0,0.97)
    save(fig,'figS7_gallery')
def figS7():
    G=[('LuPc2bare_gas','LuPc₂, 45° (optimised)'),('LuPc2bare_rot35','LuPc₂ rotated to 35°'),('LuPc2bare_rot25','LuPc₂ rotated to 25°'),('LuPc8_td','Lu-8, optimised')]
    fig=plt.figure(figsize=(7.0,2.3)); gs=GridSpec(1,4,figure=fig,wspace=0.05)
    for k,(i,t) in enumerate(G):
        ax=fig.add_subplot(gs[0,k]); off=img(ax,'gal_S7_%s.png'%i); f=os.path.join(EXI,'gal_S7_%s.json'%i)
        if os.path.exists(f) and off:
            j=json.load(open(f)); theta_arc(ax,'gal_S7_%s.png'%i,j,off,lw=1.3); ax.text(0.5,-0.02,'%s\nθ = %.1f°'%(t,j['geo']['th']),transform=ax.transAxes,ha='center',va='top',fontsize=6.6,color=INK)
        harf(ax,'abcd'[k],0.0,0.97)
    save(fig,'figS8_skew')

def _orb_read(f):
    R=[l.split() for l in open(f) if not l.startswith('#')]
    Z=np.array([int(r[0]) for r in R]); X=np.array([[float(x) for x in r[1:4]] for r in R])
    V=np.full((len(R),4),np.nan)
    for i,r in enumerate(R):
        if len(r)==8: V[i]=[float(x) for x in r[4:8]]
    return Z,X,V
OP,ON='#2f6fd6','#e0623a'   # orbital sign: positive, negative
ORB_FILE='orbitals_a1u_LuPc2opt.txt'   # beta HOMO/LUMO values on the C/N atoms (Figure 1d)
def ring_top(ax,Z,X,V,ring,col):
    """a1u lobes of one ring seen along the axis, on the face turned to the other ring (beta HOMO)"""
    from matplotlib.patches import Circle
    s=X[:,2]>0 if ring==1 else X[:,2]<0; k=1 if ring==1 else 0   # ring 1 = upper ring (blue in panels a-c)
    D=np.linalg.norm(X[:,None]-X[None],axis=2)
    for i in np.where(s)[0]:
        for j in np.where(s)[0]:
            if j>i and Z[i]>1 and Z[j]>1 and D[i,j]<1.6: ax.plot(X[[i,j],0],X[[i,j],1],color=col,lw=0.7,alpha=0.55,zorder=1)
    for i in np.where(s)[0]:
        v=V[i,k]
        if Z[i] in (6,7) and abs(v)>0.004: ax.add_patch(Circle(X[i,:2],0.95*(abs(v)/0.05)**0.5,fc=OP if v>0 else ON,ec='white',lw=0.3,alpha=0.9,zorder=2))
    ax.add_patch(Circle((0,0),0.5,fc='#2a9d6f',ec='white',lw=0.4,zorder=3)); ax.set_xlim(-8,8); ax.set_ylim(-8,8); ax.set_aspect('equal'); ax.axis('off')
def pair(ax,X,V,i1,i2,orb):
    """p lobes of an eclipsed C-alpha pair (ring 1 above, ring 2 below) from the computed values"""
    from matplotlib.patches import Ellipse,Circle
    k=0 if orb=='HOMO' else 2
    for y,i,c in ((2.6,i1,'#1d3f9e'),(0,i2,'#a8541a')):
        for s,v in ((1,V[i,k]),(-1,V[i,k+1])):
            h=0.95*(abs(v)/0.05)**0.5; ax.add_patch(Ellipse((0,y+s*h*0.55),0.62*h,1.1*h,fc=OP if v>0 else ON,ec='white',lw=0.4,alpha=0.92))
        ax.add_patch(Circle((0,y),0.11,fc=c,ec='white',lw=0.3,zorder=3))
    ax.set_xlim(-0.8,0.8); ax.set_ylim(-1.2,3.8); ax.set_aspect('equal'); ax.axis('off')
def panel_d(fig,ss):
    D=fig.add_subplot(ss); D.axis('off'); D.set_xlim(0,10); D.set_ylim(0,10); harf(D,'d',0.0,0.98)
    c1,c2,cm='#1d3f9e','#a8541a',INK
    D.plot([0.5,2.3],[5,5],color=c1,lw=2.2); D.plot([7.7,9.5],[5,5],color=c2,lw=2.2)
    D.text(1.4,5.35,'a$_{1u}$ (ring 1)',ha='center',va='bottom',fontsize=7,color=c1); D.text(8.6,5.35,'a$_{1u}$ (ring 2)',ha='center',va='bottom',fontsize=7,color=c2)
    D.plot([3.9,6.1],[8.4,8.4],color=cm,lw=2.2); D.plot([3.9,6.1],[1.6,1.6],color=cm,lw=2.2)
    for xa,xb,ya,yb in ((2.3,3.9,5,8.4),(2.3,3.9,5,1.6),(7.7,6.1,5,8.4),(7.7,6.1,5,1.6)): D.plot([xa,xb],[ya,yb],color=MUT,lw=0.7,ls=(0,(2,2)))
    D.text(5,8.8,'out-of-phase (β LUMO)',ha='center',va='bottom',fontsize=7,color=cm); D.text(5,1.2,'in-phase (β HOMO)',ha='center',va='top',fontsize=7,color=cm)
    D.annotate('',xy=(5.55,8.25),xytext=(5.55,1.75),arrowprops=dict(arrowstyle='-|>',color='#c2410c',lw=1.5,mutation_scale=9))
    D.text(5.75,5.4,'IV',fontsize=8,color='#c2410c',fontweight='bold',va='center'); D.text(5.75,4.55,r'$\nu = 2|V|$',fontsize=7.5,color='#c2410c',va='center')
    D.text(5,-0.5,r'$V(\theta)\ \propto\ \beta - \cos 4\theta$',ha='center',va='top',fontsize=7.5,color=INK)
    f=os.path.join(DATA,ORB_FILE)
    if not os.path.exists(f): return D
    Z,X,V=_orb_read(f)
    def ins(x0,y0,w,h): 
        p=D.get_position(); return fig.add_axes([p.x0+p.width*x0/10,p.y0+p.height*y0/10,p.width*w/10,p.height*h/10])
    ring_top(ins(-0.5,-0.2,3.5,4.9),Z,X,V,1,c1); ring_top(ins(7.0,-0.2,3.5,4.9),Z,X,V,2,c2)
    # eclipsed C-alpha pair nearest the +y axis
    r1=[i for i in range(len(Z)) if Z[i]==6 and X[i,2]>0 and abs(V[i,0])>0.04]; r2=[i for i in range(len(Z)) if Z[i]==6 and X[i,2]<0 and abs(V[i,0])>0.04]
    i1=min(r1,key=lambda i:np.hypot(X[i,0],X[i,1]-3)); i2=min(r2,key=lambda i:np.linalg.norm(X[i,:2]-X[i1,:2]))
    pair(ins(3.85,5.45,1.6,2.8),X,V,i1,i2,'LUMO'); pair(ins(3.85,1.75,1.6,2.8),X,V,i1,i2,'HOMO')
    return D

def fig1():
    """Figure 1: LuPc2 exported from the Explorer (rings coloured, Ln-N distances and N-Ln-N angle measured in the Explorer), descriptors, IV orbitals"""
    J=lambda k:json.load(open(os.path.join(EXI,'fig1_%s.json'%k))) if os.path.exists(os.path.join(EXI,'fig1_%s.json'%k)) else None
    fig=plt.figure(figsize=(7.0,5.3)); gs=GridSpec(2,2,figure=fig,width_ratios=[1.25,1],height_ratios=[1,0.9],wspace=0.08,hspace=0.12)
    A=fig.add_subplot(gs[0,0]); img(A,'fig1_a.png'); harf(A,'a',0.0,0.98)
    C=fig.add_subplot(gs[0,1]); oc=img(C,'fig1_c.png'); harf(C,'c',0.0,0.98); jc=J('c')
    if jc and oc:
        S=metal_fix('fig1_c.png',jc['scr'],oc); m=S[0]; t=S[1]; b=min(S[5:9],key=lambda q:np.linalg.norm((q-m)/np.linalg.norm(q-m)-(t-m)/np.linalg.norm(t-m)))
        for q,c in ((t,'#1d3f9e'),(b,'#a8541a')): C.plot([m[0],m[0]+3.7*(q[0]-m[0])],[m[1],m[1]+3.7*(q[1]-m[1])],color=c,lw=1.8,solid_capstyle='round')
        a0=np.degrees(np.arctan2(-(t[1]-m[1]),t[0]-m[0])); a1=np.degrees(np.arctan2(-(b[1]-m[1]),b[0]-m[0])); R=3.4*np.linalg.norm(t-m)
        lo,hi=sorted([a0,a1]); 
        if hi-lo>180: lo,hi=hi,lo+360
        tt=np.radians(np.linspace(lo,hi,40)); C.plot(m[0]+R*np.cos(tt),m[1]-R*np.sin(tt),color=INK,lw=1.2)
        tm=np.radians((lo+hi)/2); C.text(m[0]+1.2*R*np.cos(tm),m[1]-1.2*R*np.sin(tm),'θ = %.1f°'%jc['geo']['th'],fontsize=8,color=INK,ha='center',va='center',bbox=dict(boxstyle='round,pad=0.2',fc='white',ec='none',alpha=0.85))
    B=fig.add_subplot(gs[1,0]); ob=img(B,'fig1_b.png'); harf(B,'b',0.0,0.98); jb=J('b')
    if jb and ob:
        S=metal_fix('fig1_b.png',jb['scr'],ob); x0,x1=B.get_xlim(); yN=[S[9][1],S[10][1]]; yC=[S[11][1],S[12][1]]
        for y in yN: B.plot([x0,x1],[y,y],color=INK,lw=0.9,ls=(0,(4,2.5)))
        for y in yC: B.plot([x0,x1],[y,y],color='#6250d6',lw=0.9,ls=(0,(4,2.5)))
        B.annotate('',xy=(x0+0.04*(x1-x0),yN[0]),xytext=(x0+0.04*(x1-x0),yN[1]),arrowprops=dict(arrowstyle='<->',color=INK,lw=1.2,shrinkA=0,shrinkB=0))
        B.text(x0+0.06*(x1-x0),np.mean(yN),'$d_{\\mathrm{N4}}$ = %.3f Å'%jb['geo']['dN4'],fontsize=8,color=INK,va='center',ha='left',bbox=dict(boxstyle='round,pad=0.15',fc='white',ec='none',alpha=0.85))
        B.annotate('',xy=(x1-0.04*(x1-x0),yC[0]),xytext=(x1-0.04*(x1-x0),yC[1]),arrowprops=dict(arrowstyle='<->',color='#6250d6',lw=1.2,shrinkA=0,shrinkB=0))
        B.text(x1-0.06*(x1-x0),np.mean(yC),'$d_{\\mathrm{C\\alpha}}$ = %.3f Å'%jb['geo']['dCa'],fontsize=8,color='#6250d6',va='center',ha='right',bbox=dict(boxstyle='round,pad=0.15',fc='white',ec='none',alpha=0.85))
        B.set_xlim(x0,x1)
    panel_d(fig,gs[1,1])
    save(fig,'fig1_structure')

CG='#2a78d6'; CK='#eb6834'; CI=INK   # orbital gap, electron-hole term, IV band
def fig6():
    R={r['id']:r for r in csv.DictReader(open(os.path.join(DATA,'aryl_series.csv')))}; T=[R[k] for k in ('t35','xray','t70','t90')]; c0,me=R['core'],R['me8']
    x=np.array([float(r['cos2_mean']) for r in T]); ev=lambda r:float(r['IV_cm-1'])/8065.544; xx=np.linspace(0,0.75,50)
    fig=plt.figure(figsize=(7.0,5.0)); gs=GridSpec(2,3,figure=fig,wspace=0.48,hspace=0.36,height_ratios=[0.9,1]); A,B,C=[fig.add_subplot(gs[1,i]) for i in range(3)]
    for k,(f,t) in enumerate((('aryl_t35.png','phenyl torsion 35°'),('aryl_xray.png','X-ray torsions'),('aryl_t90.png','phenyl torsion 90°'))):
        ax=fig.add_subplot(gs[0,k]); img(ax,f,(0.2,0.8,0.0,0.82)); ilab(ax,t,-0.02); harf(ax,'abc'[k],-0.05,0.98)
    y=np.array([float(r['shift_vs_core_pct']) for r in T]); b,a=np.polyfit(x,y,1)
    A.axhline(0,color=HAIR,lw=0.7); A.plot(xx,a+b*xx,color=CI,lw=1.1); A.scatter(x,y,s=26,c=CI,edgecolors='white',linewidths=0.8,zorder=3)
    for xi,yi,r in zip(x,y,T): A.annotate('%.0f°'%float(r['torsion_mean_deg']),(xi,yi),xytext=(6,4),textcoords='offset points',fontsize=6.2,color=INK2)
    A.scatter([0],[float(me['shift_vs_core_pct'])],s=26,marker='s',facecolors='white',edgecolors=INK2,linewidths=0.9,zorder=3)
    A.annotate('methyl',(0,float(me['shift_vs_core_pct'])),xytext=(6,0),textcoords='offset points',fontsize=6.2,color=INK2,va='center')
    A.set_xlim(-0.04,0.75); A.set_ylim(-5.8,0.4); A.set_xlabel('mean cos²φ of the phenyl torsion'); A.set_ylabel('shift of the IV band vs core (%)'); ygrid(A); harf(A,'d',-0.3,1.02)
    D={'IV band':(np.array([ev(r)-ev(c0) for r in T]),CI,'o'),'orbital gap':(np.array([float(r['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV']) for r in T]),CG,'s'),
       'electron–hole term $K$':(np.array([float(r['K_eV'])-float(c0['K_eV']) for r in T]),CK,'D')}
    mg=float(me['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV']); mv=ev(me)-ev(c0); MEv={'IV band':mv,'orbital gap':mg,'electron–hole term $K$':mv-mg}
    B.axhline(0,color=HAIR,lw=0.7)
    for lab,(yy,c,mk) in D.items():
        b,a=np.polyfit(x,1000*yy,1); B.plot(xx,a+b*xx,color=c,lw=1.1); B.scatter(x,1000*yy,s=20,marker=mk,c=c,edgecolors='white',linewidths=0.7,zorder=3,label=lab)
        B.scatter([{'IV band':-0.12,'orbital gap':-0.1}.get(lab,-0.08)],[1000*MEv[lab]],s=20,marker=mk,facecolors='white',edgecolors=c,linewidths=0.9,zorder=3)
    B.axvline(-0.05,color=HAIR,lw=0.7); B.text(-0.135,-30,'methyl',fontsize=6.0,color=INK2,ha='left',va='top')
    B.set_xlim(-0.15,0.75); B.set_xticks([0,0.2,0.4,0.6]); B.set_xlabel('mean cos²φ'); B.set_ylabel('change vs core (meV)'); B.legend(loc='lower left',bbox_to_anchor=(0.0,0.0),fontsize=6.0,handletextpad=0.3); ygrid(B); harf(B,'e',-0.3,1.02)
    for k,c,mk,lab in (('aryl_share_beta_HOMO',CI,'o','β HOMO'),('aryl_share_beta_LUMO',MUT,'o','β LUMO')):
        yy=np.array([float(r[k]) for r in T]); b,a=np.polyfit(x,yy,1); C.plot(xx,a+b*xx,color=c,lw=1.0); C.scatter(x,yy,s=20,c=c,edgecolors='white',linewidths=0.7,zorder=3,label=lab)
    lo,hi=sorted([float(me['aryl_share_beta_LUMO']),float(me['aryl_share_beta_HOMO'])]); C.axhspan(lo,hi,color=WASH,lw=0,zorder=0)
    C.text(0.74,lo-0.0015,'methyl at the same positions',ha='right',va='top',fontsize=6.0,color=INK2)
    C.set_xlim(-0.04,0.75); C.set_ylim(0,0.06); C.set_xlabel('mean cos²φ'); C.set_ylabel('phenyl share of the orbital'); C.legend(loc='upper left',fontsize=6.2); ygrid(C); harf(C,'f',-0.3,1.02)
    save(fig,'fig6_aryl')
CLS=['N-linked cyclic diarylamine','thioether','C-linked heteroaryl','alkyl','alkoxy','amino']
def fig7():
    R=list(csv.DictReader(open(os.path.join(DATA,'heteroleptic.csv')))); cC=json.load(open(os.path.join(DATA,'constants.json')))['Ca_form']['donor_c']
    fig=plt.figure(figsize=(7.0,5.2)); gs=GridSpec(2,2,figure=fig,width_ratios=[1.15,1],wspace=0.32,hspace=0.42,height_ratios=[0.95,1]); A,B=fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])
    SP={r['id']:r for r in R}
    for k,(f,id,t) in enumerate((('het_ptz.png','LuHt4_td2','phenothiazin-10-yl (N-linked)'),('het_nme2.png','LuHtNMe2_opt_td','dimethylamino'))):
        ax=fig.add_subplot(gs[0,k]); img(ax,f); hd=float(SP[id]['hole_on_donor_ring']); ilab(ax,'%s\nhole: substituted ring (top) %.2f, unsubstituted ring %.2f'%(t,hd,1-hd),-0.02); harf(ax,'ab'[k],-0.02,0.98)
    A.axvline(0.5,color=INK2,lw=0.6,ls=(0,(2.5,2)))
    for i,c in enumerate(CLS):
        rr=[r for r in R if r['donor_class']==c]
        for j,r in enumerate(rr):
            A.scatter(float(r['hole_on_donor_ring']),i+0.12*(j-(len(rr)-1)/2),s=26,marker=MK[r['metal']],c=MC[r['metal']],edgecolors='white',linewidths=0.7,zorder=3)
    A.set_yticks(range(len(CLS))); A.set_yticklabels(CLS); A.set_ylim(-0.6,len(CLS)-0.4); A.invert_yaxis(); A.set_xlim(0.33,0.72)
    A.set_xlabel('share of the ground-state hole on the substituted ring'); A.grid(axis='x',color=HAIR,lw=0.45); A.set_axisbelow(True)
    A.text(0.49,len(CLS)-0.45,'hole on the\nunsubstituted ring',ha='right',va='bottom',fontsize=6.0,color=INK2); A.text(0.51,len(CLS)-0.45,'hole on the\nsubstituted ring',ha='left',va='bottom',fontsize=6.0,color=INK2)
    for m in ('Lu','Y','La'): A.scatter([],[],marker=MK[m],c=MC[m],s=22,label=m)
    A.legend(loc='upper right',ncol=3,fontsize=6.2,bbox_to_anchor=(1.0,1.1)); harf(A,'c',-0.62,1.02)
    aa=np.linspace(0,0.4,60); B.plot(aa,100*(1/np.sqrt(1-aa**2)-1),color=INK2,lw=1.0,ls=(0,(3,2))); B.axhline(0,color=HAIR,lw=0.7)
    B.text(0.385,100*(1/np.sqrt(1-0.385**2)-1)+0.5,'1/√(1 − $a^2$)',ha='right',fontsize=6.2,color=INK2)
    for r in R:
        a=abs(float(r['ring_spin_asym_a'])); dv=100*((1+float(r['dev_Ca_pct'])/100)/(1+cC*float(r['D_sub']))-1)
        B.scatter(a,dv,s=26,marker=MK[r['metal']],c=MC[r['metal']],edgecolors='white',linewidths=0.7,zorder=3)
        if r['donor'] in ('dimethylamino','octakis(methylthio)') or r['metal']=='La': B.annotate({'octakis(methylthio)':'Lu, octakis-\n(methylthio)'}.get(r['donor'],r['metal']+', '+r['donor']),(a,dv),xytext={'dimethylamino':(-4,-9),'octakis(methylthio)':(3,6)}.get(r['donor'],(5,5)),textcoords='offset points',ha='center' if r['donor']=='dimethylamino' else 'left',va='top' if r['donor']=='dimethylamino' else 'bottom',fontsize=5.8,color=INK2)
    B.set_xlim(0,0.4); B.set_xlabel('ring spin asymmetry |$a$|'); B.set_ylabel('deviation from the donor-corrected\nCα form (%)'); ygrid(B); harf(B,'d',-0.25,1.02)
    save(fig,'fig7_heteroleptic')
def capture(fn,name):
    b=io.StringIO()
    with contextlib.redirect_stdout(b): fn()
    open(os.path.join(OUT,name),'w').write(b.getvalue()); print('wrote',name)
if __name__=='__main__':
    fig1(); P=table1(); fig2(P); fig3(P); fig6(); fig7(); figS1(P); figS2(); figS3(); figS4(); figS5(); figS6(); figS7()
    import aryl,heteroleptic
    capture(aryl.series,'tableS_aryl.txt'); capture(heteroleptic.table,'tableS_heteroleptic.txt')
    import subprocess; open(os.path.join(OUT,'tableS_unsubstituted_LnPc2.txt'),'w').write(subprocess.run([sys.executable,os.path.join(H,'unsubstituted_LnPc2_test.py')],capture_output=True,text=True).stdout); print('wrote tableS_unsubstituted_LnPc2.txt')
