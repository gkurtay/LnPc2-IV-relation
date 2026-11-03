"""Draw Figures 1-3, 6 and 7 and Figure S4 in the layout of the published article (white panels, no image cards),
from the files in data/ and the Explorer images in outputs/explorer/.

usage: python article_figures.py [output_directory]          (default: ../outputs/article)
needs: numpy, scipy, matplotlib (as make_outputs.py)

make_outputs.py reproduces the same tables and the same data in its own layout; this script only redraws the figures
as they appear in the article. Figures 4 and 5 combine Explorer images with literature spectra and are not regenerated.
Panel letters only; captions are in the article."""
import sys,os,io,contextlib
H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,H)
OUTA=os.path.abspath(sys.argv[1] if len(sys.argv)>1 else os.path.join(H,'..','outputs','article')); os.makedirs(OUTA,exist_ok=True)
sys.argv=['article_figures.py',OUTA]
with contextlib.redirect_stdout(io.StringIO()):
    import make_outputs as M
    from make_outputs import *
    P=M.table1()
M.EXI=os.path.join(H,'..','outputs','explorer'); EXI=M.EXI   # structure panels exported from the Explorer, as deposited

# ---- shared drawing helpers (halo text, dumbbell bars, statistics line, bands)
"""Ortak tasarim yardimcilari (v30): hale, dumbbell cubugu, istatistik satiri, dogrudan etiket."""
import numpy as np, matplotlib.patheffects as pe
BGc='white'; INKc='#1d1d1b'; INK2c='#55554f'; GRc='#8a94a6'
def halo(t,w=2.6,c=BGc): t.set_path_effects([pe.withStroke(linewidth=w,foreground=c)]); return t
def htext(ax,x,y,s,w=2.6,**k): return halo(ax.text(x,y,s,**k),w)
def bar(ax,x,y0,y1,col,lw=4.2,alpha=0.18,z=2): ax.vlines([x],[y0],[y1],color=col,lw=lw,alpha=alpha,zorder=z,capstyle='round')
def stat(ax,s,x=0.97,y=0.96,ha='right',va='top',size=6.2,col=INK2c):
    return halo(ax.text(x,y,s.replace('-','−'),transform=ax.transAxes,ha=ha,va=va,fontsize=size,color=col,linespacing=1.3,zorder=9))
def band(ax,y0,y1,col=GRc,alpha=0.2,x=None,z=1):
    if x is None: return ax.axhspan(y0,y1,color=col,alpha=alpha,lw=0,zorder=z)
    return ax.fill_between(x,y0,y1,color=col,alpha=alpha,lw=0,zorder=z)

import numpy as np, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
NEED=[]

# ---- Figure 1
M.fig1()

# ---- Figure 2
from matplotlib.gridspec import GridSpecFromSubplotSpec as SGS
p=P['Ca']; node=np.degrees(np.arccos(p[2]))/4; DC=fl('LuPc2bare_gas','d_Ca_A'); rng=np.random.default_rng(7); tt=np.linspace(node,45,160); B=[]
for _ in range(1000):
    k=rng.choice(len(nu),len(nu)); q=fitq(gC,dc[k],t0[k],nu[k],p); B.append(gC(q,DC,tt))
lo,hi=np.percentile(np.array(B),[2.5,97.5],axis=0)
R=sorted([k for k in S if k.startswith('LuPc2bare_rot') and k!='LuPc2bare_rot45']+['LuPc2bare_gas'],key=lambda k:fl(k,'theta_deg'))
th=np.array([fl(k,'theta_deg') for k in R]); v=np.array([fl(k,'IV_cm-1') for k in R]); f=np.array([fl(k,'f_osc') for k in R])
a=np.abs(np.array([float(S[k]['ring_spin_asym_a'] or 'nan') for k in R]))
fig=plt.figure(figsize=(7.0,5.0)); gs=GridSpec(3,3,figure=fig,width_ratios=[1,1,0.9],height_ratios=[0.95,1,1],wspace=0.38,hspace=0.42)
g0=SGS(1,3,subplot_spec=gs[0,0:2],wspace=0.02)
for k,(t,lab) in enumerate((('0','eclipsed'),('20','coupling node'),('45','staggered'))):
    ax=fig.add_subplot(g0[0,k]); off=img(ax,'rot_top_%s.png'%t); jf=os.path.join(EXI,'rot_top_%s.json'%t)
    if os.path.exists(jf) and off and t!='0': theta_arc(ax,'rot_top_%s.png'%t,json.load(open(jf)),off,lw=1.2)
    th_=json.load(open(jf))['geo']['th'] if os.path.exists(jf) else float(t)
    ax.text(0.5,-0.03,'θ = %.0f°, %s'%(abs(th_),lab),transform=ax.transAxes,ha='center',va='top',fontsize=6.6,color=INK)
    if k==0: harf(ax,'a',-0.05,0.98)
# (b)
A=fig.add_subplot(gs[1:,0:2]); ygrid(A)
GRc_='#8a94a6'
A.axvspan(-1.5,node,color=GRc_,alpha=0.10,lw=0,zorder=0); A.axvspan(24.5,46.5,color=MC['Lu'],alpha=0.07,lw=0,zorder=0)
A.axvline(node,color=INK2,lw=0.7,ls=(0,(2.5,2)),zorder=1)
A.text(node,11780,'node %.1f°'%node,ha='center',va='center',fontsize=6.4,color='white',fontweight='bold',bbox=dict(boxstyle='round,pad=0.3,rounding_size=0.5',fc=INK2,ec='none'),zorder=8)
htext(A,9.2,11000,'below the node',ha='center',fontsize=6.4,color=INK2); htext(A,35.5,11000,'fitted range',ha='center',fontsize=6.4,color=MC['Lu'])
A.fill_between(tt,lo,hi,color=MC['Lu'],alpha=0.16,lw=0,zorder=2); A.plot(tt,gC(p,DC,tt),color=MC['Lu'],lw=1.6,zorder=3)
m=th>=24; A.scatter(th[m],v[m],s=30,c=MC['Lu'],edgecolors='white',linewidths=0.8,zorder=4)
lo_=(th<24)&(np.abs(th-20)>0.5); A.scatter(th[lo_],v[lo_],s=28,facecolors='white',edgecolors=INK2,linewidths=0.9,zorder=4)
k20=int(np.argmin(np.abs(th-20))); A.scatter([th[k20]],[v[k20]],s=32,marker='D',c=INK2,edgecolors='white',linewidths=0.8,zorder=4)
halo(A.annotate('ring-localised ground state\n|$a$| = %.2f, $f$ = %.3f'%(a[k20],f[k20]),xy=(th[k20],v[k20]),xytext=(1.0,700),fontsize=6.3,color=INK,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.6,shrinkA=0,shrinkB=3)))
i_=95; ang_=np.degrees(np.arctan2((gC(p,DC,tt[i_+3])-gC(p,DC,tt[i_-3]))/12500*A.bbox.height,(tt[i_+3]-tt[i_-3])/48*A.bbox.width))
htext(A,tt[i_]+1.2,gC(p,DC,tt[i_])-650,'Cα form, 95% band',fontsize=6.3,color=MC['Lu'],rotation=ang_,rotation_mode='anchor',ha='center',va='top')
A.set_xlim(-1.5,46.5); A.set_xticks(range(0,46,5)); A.set_xlabel('skew angle θ (°)'); A.set_ylabel('IV energy (cm$^{-1}$)'); A.set_ylim(0,12400); A.set_yticks(range(0,12001,2000)); harf(A,'b',-0.1,1.02)
# (c)
Bf=fig.add_subplot(gs[0,2]); ygrid(Bf); Bf.axvline(node,color=INK2,lw=0.7,ls=(0,(2.5,2)))
Bf.fill_between(th,0,f,color=MC['Lu'],alpha=0.10,lw=0,zorder=1); Bf.plot(th,f,color=MC['Lu'],lw=1.2,zorder=2); Bf.scatter(th,f,s=13,c=MC['Lu'],edgecolors='white',linewidths=0.5,zorder=3)
halo(Bf.annotate('$f$ = %.3f'%f[k20],(th[k20],f[k20]),xytext=(0,30),textcoords='offset points',fontsize=6,color=INK2,ha='center',arrowprops=dict(arrowstyle='-',color=INK2,lw=0.5,shrinkA=0,shrinkB=2)))
Bf.set_xlim(-1.5,46.5); Bf.set_xticks(range(0,46,15)); Bf.set_ylim(0,0.18); Bf.set_ylabel('oscillator strength $f$'); Bf.set_xlabel('θ (°)',labelpad=1); harf(Bf,'c',-0.36,1.04)
g1=SGS(2,1,subplot_spec=gs[1:,2],hspace=0.25)
for k,(t,kk) in enumerate((('20',k20),('25',int(np.argmin(np.abs(th-25)))))):
    ax=fig.add_subplot(g1[k,0]); img(ax,'rot_spin_%s.png'%t)
    ax.text(0.5,-0.02,'θ = %.1f°, |$a$| = %.2f'%(th[kk],a[kk]),transform=ax.transAxes,ha='center',va='top',fontsize=6.6,color=INK)
    if k==0: harf(ax,'d',-0.12,0.98)
save(fig,'fig2_rotation_scan')

# ---- Figure 3
fig=plt.figure(figsize=(7.0,3.8)); gs=GridSpec(2,2,figure=fig,height_ratios=[3.2,1.05],hspace=0.08,wspace=0.3); met=[S[k]['metal'] for k in FIT]
rot=np.array(['rot' in k for k in FIT])
for j,(fm,lab,s) in enumerate((('Ca','Cα form','a'),('N4','N$_4$ form','b'))):
    g,d,_=FORMS[fm]; L=np.array([g(fitq(g,d[np.arange(len(nu))!=k],t0[np.arange(len(nu))!=k],nu[np.arange(len(nu))!=k],P[fm]),d[k],t0[k]) for k in range(len(nu))])
    a=fig.add_subplot(gs[0,j]); r=fig.add_subplot(gs[1,j],sharex=a); lim=(2400,10900); e=np.sqrt(np.mean((L-nu)**2))
    def draw(ax,big=True):
        ax.fill_between(lim,[lim[0]-e,lim[1]-e],[lim[0]+e,lim[1]+e],color=MC['Lu'],alpha=0.12,lw=0,zorder=0); ax.plot(lim,lim,color=INK2,lw=0.6,ls=(0,(2.5,2)),zorder=1)
        for k in range(len(nu)):
            m=met[k]; kw=dict(s=36 if big else 30,marker=MK[m],zorder=3)
            if rot[k]: ax.scatter([nu[k]],[L[k]],facecolors='white',edgecolors=MC[m],linewidths=1.0,**kw)
            else: ax.scatter([nu[k]],[L[k]],c=MC[m],edgecolors='white',linewidths=0.7,**kw)
    draw(a)
    for k in range(len(nu)):
        if rot[k] and nu[k]<9000: htext(a,nu[k]+170,L[k]-120,SHORT[FIT[k]],fontsize=6,color=INK2,va='top')
        if FIT[k]=='LaPc2bare_td': htext(a,nu[k]+200,L[k]-120,'LaPc$_2$, La-1',fontsize=6,color=INK2,ha='left',va='top')
    a.set_xlim(lim); a.set_ylim(lim); a.set_ylabel('LOO prediction (cm$^{-1}$)'); plt.setp(a.get_xticklabels(),visible=False); a.tick_params(axis='x',length=0)
    chip(a,lab,MC['Lu'],x=0.97,y=0.30,size=7.0); a.texts[-1].set_ha('right')
    a.text(0.97,0.17,'%.0f'%e,transform=a.transAxes,ha='right',va='top',fontsize=15,fontweight='bold',color=INK)
    a.text(0.97,0.055,'cm$^{-1}$ LOO RMS · n = %d'%len(nu),transform=a.transAxes,ha='right',va='top',fontsize=6.2,color=INK2)
    ia=a.inset_axes([0.03,0.54,0.42,0.44]); ia.set_facecolor('white'); ia.white=True; draw(ia,False); c0,c1=9230,10400; ia.set_xlim(c0,c1); ia.set_ylim(c0,c1)
    ia.set_xticks([]); ia.set_yticks([])
    for sp in ia.spines.values(): sp.set_visible(True); sp.set_color(MUT); sp.set_linewidth(0.5)
    ks=[k for k in range(len(nu)) if c0<nu[k]<c1]
    for side in (1,-1):
        kk=sorted([k for k in ks if (1 if L[k]<=nu[k] else -1)==side],key=lambda k:L[k]); ys=[L[k] for k in kk]
        for it in range(50):
            for q in range(1,len(ys)):
                if ys[q]-ys[q-1]<62: d_=(62-(ys[q]-ys[q-1]))/2; ys[q-1]-=d_; ys[q]+=d_
        for k,y in zip(kk,ys):
            ia.annotate(SHORT[FIT[k]],(nu[k],L[k]),xytext=(nu[k]+side*(70 if side>0 else 45),y),textcoords='data',fontsize=5.2,color=INK2,ha='left' if side>0 else 'right',va='center')
    a.indicate_inset_zoom(ia,edgecolor=MUT,alpha=0.6,lw=0.5)
    # residuals as stems
    r.axhspan(-e,e,color=MC['Lu'],alpha=0.12,lw=0); r.axhline(0,color=INK2,lw=0.6)
    for k in range(len(nu)):
        m=met[k]; res=L[k]-nu[k]; r.plot([nu[k],nu[k]],[0,res],color=MC[m],lw=1.0,alpha=0.7,zorder=2,solid_capstyle='round')
        kw=dict(s=22,marker=MK[m],zorder=3)
        if rot[k]: r.scatter([nu[k]],[res],facecolors='white',edgecolors=MC[m],linewidths=1.0,**kw)
        else: r.scatter([nu[k]],[res],c=MC[m],edgecolors='white',linewidths=0.7,**kw)
    htext(r,2550,e+25,'± LOO RMS',fontsize=5.6,color=MC['Lu'],ha='left',va='bottom')
    r.set_ylim(-520,520); r.set_yticks([-400,0,400]); r.set_ylabel('Δ (cm$^{-1}$)'); r.set_xlabel('computed IV energy (cm$^{-1}$)'); a.spines['bottom'].set_visible(False); ygrid(r); harf(a,s,-0.17,1.02)
H=[Line2D([],[],marker=MK[m],ls='',mfc=MC[m],mec='white',ms=5.5,label=m+'Pc$_2$ family, optimised') for m in ('Lu','Y','La')]+[Line2D([],[],marker='o',ls='',mfc='white',mec=MC['Lu'],mew=1.0,ms=5,label='LuPc$_2$, rigid rotation')]
fig.legend(handles=H,loc='upper center',ncol=4,bbox_to_anchor=(0.52,1.03),fontsize=6.3,columnspacing=1.2,handletextpad=0.3)
save(fig,'fig3_loo')

# ---- Figure 6
R6={r['id']:r for r in csv.DictReader(open(os.path.join(DATA,'aryl_series.csv')))}; T6=[R6[k] for k in ('t35','xray','t70','t90')]; c0,me=R6['core'],R6['me8']
x=np.array([float(r['cos2_mean']) for r in T6]); ev=lambda r:float(r['IV_cm-1'])/8065.544; xx=np.linspace(-0.02,0.74,50)
tors=[float(r['torsion_mean_deg']) for r in T6]
fig=plt.figure(figsize=(7.0,5.1)); gs=GridSpec(2,3,figure=fig,wspace=0.5,hspace=0.42,height_ratios=[0.9,1]); A,B,C=[fig.add_subplot(gs[1,i]) for i in range(3)]
for k,(f,t) in enumerate((('aryl_t35.png','phenyl torsion 35°'),('aryl_xray.png','X-ray torsions'),('aryl_t90.png','phenyl torsion 90°'))):
    ax=fig.add_subplot(gs[0,k]); img(ax,f,(0.2,0.8,0.0,0.82)); ilab(ax,t,-0.02); harf(ax,'abc'[k],-0.05,0.98)
def top(ax):
    sx=ax.secondary_xaxis('top'); sx.set_xticks(x); sx.set_xticklabels(['%.0f°'%t for t in tors],fontsize=5.8,color=INK2); sx.tick_params(length=0,pad=1.5)
    for sp in sx.spines.values(): sp.set_visible(False)
    return sx
METc='#8a94a6'
# (d)
y=np.array([float(r['shift_vs_core_pct']) for r in T6]); b,a=np.polyfit(x,y,1); ym=float(me['shift_vs_core_pct'])
ygrid(A); A.axhline(0,color=INK2,lw=0.6,ls=(0,(2.5,2))); A.axhspan(ym-0.07,ym+0.07,color=METc,alpha=0.35,lw=0,zorder=1)
htext(A,0.73,ym+0.12,'methyl',fontsize=5.9,color=INK2,ha='right',va='bottom')
A.plot(xx,a+b*xx,color=CI,lw=1.2,zorder=2); A.scatter(x,y,s=30,c=CI,edgecolors='white',linewidths=0.8,zorder=3)
stat(A,'slope %.1f%% per unit\ncos²φ'%b,x=0.04,y=0.05,ha='left',va='bottom')
A.set_xlim(-0.04,0.75); A.set_ylim(-5.8,0.4); A.set_xlabel('mean cos²φ'); A.set_ylabel('shift of the IV band vs core (%)'); top(A).set_xlabel('mean phenyl torsion φ',fontsize=6.2,color=INK2,labelpad=3); harf(A,'d',-0.3,1.13)
# (e)
D={'IV band':(np.array([ev(r)-ev(c0) for r in T6]),CI,'o'),'orbital gap':(np.array([float(r['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV']) for r in T6]),CG,'s'),
   'electron–hole\nterm $K$':(np.array([float(r['K_eV'])-float(c0['K_eV']) for r in T6]),CK,'D')}
mg=float(me['beta_HOMO_LUMO_gap_eV'])-float(c0['beta_HOMO_LUMO_gap_eV']); mv=ev(me)-ev(c0); MEv={'IV band':mv,'orbital gap':mg,'electron–hole\nterm $K$':mv-mg}
ygrid(B); B.axhline(0,color=INK2,lw=0.6,ls=(0,(2.5,2)))
B.axvspan(-0.2,-0.045,color=METc,alpha=0.16,lw=0,zorder=0); htext(B,-0.1225,-75,'methyl',fontsize=5.9,color=INK2,ha='center',va='top',rotation=90)
END={'IV band':(0,-2),'orbital gap':(0,-2),'electron–hole\nterm $K$':(0,4)}
for lab,(yy,c,mk) in D.items():
    b,a=np.polyfit(x,1000*yy,1); B.plot(xx[xx>=0],a+b*xx[xx>=0],color=c,lw=1.2,zorder=2); B.scatter(x,1000*yy,s=22,marker=mk,c=c,edgecolors='white',linewidths=0.7,zorder=3)
    B.scatter([{'IV band':-0.17,'orbital gap':-0.12}.get(lab,-0.07)],[1000*MEv[lab]],s=20,marker=mk,facecolors='white',edgecolors=c,linewidths=0.9,zorder=3)
    lx,ly,va_={'IV band':(0.5,a+b*0.5+7,'bottom'),'orbital gap':(0.40,a+b*0.40-8,'top'),'electron–hole\nterm $K$':(0.0,24,'center')}[lab]
    htext(B,lx,ly,lab.replace('\n',' '),fontsize=6.0,color=c,ha='left',va=va_)
B.set_xlim(-0.2,0.75); B.set_xticks([0,0.2,0.4,0.6]); B.set_xlabel('mean cos²φ'); B.set_ylabel('change vs core (meV)'); top(B); harf(B,'e',-0.3,1.13)
# (f)
ygrid(C)
for k,c,lab in (('aryl_share_beta_HOMO',CI,'β HOMO'),('aryl_share_beta_LUMO',MUT,'β LUMO')):
    yy=np.array([float(r[k]) for r in T6]); b,a=np.polyfit(x,yy,1); C.plot(xx,a+b*xx,color=c,lw=1.2); C.scatter(x,yy,s=22,c=c,edgecolors='white',linewidths=0.7,zorder=3)
    htext(C,0.76,a+b*0.74,lab,fontsize=6.0,color=c if c!=MUT else INK2,ha='left',va='center')
lo,hi=sorted([float(me['aryl_share_beta_LUMO']),float(me['aryl_share_beta_HOMO'])]); C.axhspan(lo,hi,color=METc,alpha=0.35,lw=0,zorder=1)
htext(C,0.73,lo-0.0015,'methyl at the same positions',ha='right',va='top',fontsize=5.9,color=INK2)
C.set_xlim(-0.04,0.75); C.set_ylim(0,0.06); C.set_xlabel('mean cos²φ'); C.set_ylabel('phenyl share of the orbital'); top(C); harf(C,'f',-0.3,1.13)
save(fig,'fig6_aryl')

# ---- Figure 7
R7=list(csv.DictReader(open(os.path.join(DATA,'heteroleptic.csv')))); cC=json.load(open(os.path.join(DATA,'constants.json')))['Ca_form']['donor_c']
fig=plt.figure(figsize=(7.0,5.2)); gs=GridSpec(2,2,figure=fig,width_ratios=[1.15,1],wspace=0.32,hspace=0.42,height_ratios=[0.95,1]); A,B=fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])
SP={r['id']:r for r in R7}
for k,(f,id,t) in enumerate((('het_ptz.png','LuHt4_td2','phenothiazin-10-yl (N-linked)'),('het_nme2.png','LuHtNMe2_opt_td','dimethylamino'))):
    ax=fig.add_subplot(gs[0,k]); img(ax,f); hd=float(SP[id]['hole_on_donor_ring']); ilab(ax,'%s\nhole: substituted ring (top) %.2f, unsubstituted ring %.2f'%(t,hd,1-hd),-0.02); harf(ax,'ab'[k],-0.02,0.98)
# (c) zebra satirlar, iki yari
n=len(CLS)
A.axvspan(0.33,0.5,color='#a8541a',alpha=0.04,lw=0,zorder=0); A.axvspan(0.5,0.72,color=MC['Lu'],alpha=0.04,lw=0,zorder=0)
A.axvline(0.5,color=INK2,lw=0.7,ls=(0,(2.5,2)))
for i,c in enumerate(CLS):
    rr=sorted([r for r in R7 if r['donor_class']==c],key=lambda r:float(r['hole_on_donor_ring']))
    xs_=[float(r['hole_on_donor_ring']) for r in rr]
    for j,r in enumerate(rr):
        A.scatter(float(r['hole_on_donor_ring']),i+0.13*(j-(len(rr)-1)/2),s=28,marker=MK[r['metal']],c=MC[r['metal']],edgecolors='white',linewidths=0.7,zorder=3)
A.set_yticks(range(n)); A.set_yticklabels(CLS); A.set_ylim(-0.5,n-0.5); A.invert_yaxis(); A.set_xlim(0.33,0.72)
A.set_xlabel('share of the ground-state hole on the substituted ring'); A.grid(axis='x',color='#e8eaef',lw=0.55); A.grid(axis='y',color='#e8eaef',lw=0.55); A.set_axisbelow(True)
A.text(0.49,-0.75,'← unsubstituted ring',ha='right',va='bottom',fontsize=6.0,color='#a8541a'); A.text(0.51,-0.75,'substituted ring →',ha='left',va='bottom',fontsize=6.0,color=MC['Lu'])
A.keep_style=True; A.set_facecolor(BG)
A.spines['bottom'].set_color('#b9bdc7'); A.spines['bottom'].set_linewidth(0.6); A.spines['left'].set_visible(False)
A.tick_params(axis='x',length=2.2,width=0.6,color='#b9bdc7'); A.tick_params(axis='y',length=0)
H=[Line2D([],[],marker=MK[m],ls='',mfc=MC[m],mec='white',ms=5.2,label=m) for m in ('Lu','Y','La')]
A.legend(handles=H,loc='lower right',ncol=3,fontsize=6.2,bbox_to_anchor=(1.0,1.06),handletextpad=0.2,columnspacing=0.9); harf(A,'c',-0.62,1.06)
# (d)
aa=np.linspace(0,0.4,60); cur=100*(1/np.sqrt(1-aa**2)-1)
ygrid(B); B.fill_between(aa,0,cur,color=INK2,alpha=0.07,lw=0,zorder=0); B.plot(aa,cur,color=INK2,lw=1.1,ls=(0,(3,2)),zorder=1); B.axhline(0,color=INK2,lw=0.6)
htext(B,0.392,100*(1/np.sqrt(1-0.392**2)-1)-0.25,'1/√(1 − $a^2$)',ha='right',va='top',fontsize=6.4,color=INK)
for r in R7:
    a_=abs(float(r['ring_spin_asym_a'])); dv=100*((1+float(r['dev_Ca_pct'])/100)/(1+cC*float(r['D_sub']))-1)
    B.scatter(a_,dv,s=30,marker=MK[r['metal']],c=MC[r['metal']],edgecolors='white',linewidths=0.7,zorder=3)
    if r['donor'] in ('dimethylamino','octakis(methylthio)') or r['metal']=='La':
        t_={'octakis(methylthio)':'Lu, octakis(methylthio)','dimethylamino':'Lu, dimethylamino'}.get(r['donor'],'La, '+r['donor'])
        off={'dimethylamino':(-6,-10,'right','top'),'octakis(methylthio)':(2,20,'left','bottom')}.get(r['donor'],(12,4,'left','bottom'))
        halo(B.annotate(t_,(a_,dv),xytext=off[:2],textcoords='offset points',ha=off[2],va=off[3],fontsize=5.9,color=INK2,arrowprops=dict(arrowstyle='-',color=INK2,lw=0.5,shrinkA=0,shrinkB=3)))
B.set_xlim(0,0.4); B.set_ylim(-1.2,9.6); B.set_xlabel('ring spin asymmetry |$a$|'); B.set_ylabel('deviation from the donor-corrected\nCα form (%)'); harf(B,'d',-0.25,1.06)
save(fig,'fig7_heteroleptic')

# ---- Figure S4 (symmetry-labelled orbital scan; data/orbital_scan.csv)
import glob,re,csv,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import csv as _csv
_R=list(_csv.DictReader(open(os.path.join(DATA,'orbital_scan.csv'))))
TH=[float(r['theta_deg']) for r in _R]; A=np.array([float(r['dE_A_minus_B1_anion_eV']) for r in _R]); N=np.array([float(r['dE_A_minus_B1_neutral_alpha_eV']) for r in _R])
def cross(y):
    k=[j for j in range(len(TH)-1) if y[j]>0>=y[j+1]][0]; return TH[k]+(TH[k+1]-TH[k])*y[k]/(y[k]-y[k+1])
ca,cn=cross(A),cross(N); print('kesisim anyon %.2f notr %.2f'%(ca,cn))
c4=np.cos(np.radians(4*np.array(TH))); fa=np.polyfit(c4,A,1); fn=np.polyfit(c4,N,1); print('egim',fa,fn,'r',np.corrcoef(c4,A)[0,1],np.corrcoef(c4,N)[0,1])
from figstyle import polish, MC as _MC
from matplotlib import font_manager
FONT='Inter' if any('Inter' in f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family':FONT,'font.size':7,'axes.labelsize':7.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
fig,ax=plt.subplots(1,2,figsize=(6.4,2.5)); C1,C2=_MC['Lu'],_MC['Y']
a=ax[0]; a.axhline(0,color='#999',lw=0.7); a.plot(TH,A,'o-',color=C1,ms=4,lw=1.2,label='anion (closed shell)'); a.plot(TH,N,'s--',color=C2,ms=4,mfc='white',lw=1.2,label='neutral radical (α)')
a.axvline(20.0,color='#555',ls=':',lw=0.9); a.axvline((ca+cn)/2,color='k',lw=0.9)
a.text(19.4,0.95,'IV node\n20.0°',ha='right',va='center',fontsize=7.5,color='#444'); a.text(21.8,-1.05,'orbital crossing\n%.1f°'%((ca+cn)/2),ha='left',va='center',fontsize=7.5)
a.set_xlim(0,45); a.set_xlabel('skew angle θ (°)'); a.set_ylabel('ε(A) − ε(B1) (eV)'); a.legend(frameon=False,fontsize=7.5,loc='upper right')
a.text(-0.16,1.02,'a',transform=a.transAxes,fontweight='bold',fontsize=10)
b=ax[1]; b.axhline(0,color='#999',lw=0.7); b.axvline(0,color='#999',lw=0.7)
xx=np.linspace(-1,1,10); b.plot(c4,A,'o',color=C1,ms=4); b.plot(c4,N,'s',color=C2,ms=4,mfc='white'); b.plot(xx,np.polyval(fa,xx),color=C1,lw=0.9); b.plot(xx,np.polyval(fn,xx),color=C2,lw=0.9,ls='--')
b.set_xlabel('cos 4θ'); b.set_ylabel('ε(A) − ε(B1) (eV)'); b.text(-0.16,1.02,'b',transform=b.transAxes,fontweight='bold',fontsize=10)
for q in ax: q.spines['top'].set_visible(False); q.spines['right'].set_visible(False)
fig.tight_layout(); polish(fig); fig.savefig(os.path.join(OUTA,'figS4_orbital_scan.png'),dpi=300); fig.savefig(os.path.join(OUTA,'figS4_orbital_scan.pdf')); print('ok')
