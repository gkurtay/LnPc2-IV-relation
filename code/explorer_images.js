// Structure images of the article from the LnPc2 Explorer, driven through its scripting interface (window.LNPC).
// Usage (from the package root): node code/explorer_images.js [code/explorer_jobs.json]
// Requires Node.js and Playwright (npm i playwright; npx playwright install chromium). Writes outputs/explorer/*.png
// and, next to each image, a .json with the descriptors, the measured distances/angles and the screen positions
// of the metal, the N atoms and the ring-plane centres used by make_outputs.py to place labels.
const {chromium}=require('playwright');const fs=require('fs');const path=require('path');
(async()=>{const JJ=[].concat(JSON.parse(fs.readFileSync(process.argv[2]||path.join(__dirname,'explorer_jobs.json'))));fs.mkdirSync('outputs/explorer',{recursive:true});
const b=await chromium.launch({args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
for(const J of JJ){const p=await b.newPage({viewport:{width:J.vw||1500,height:J.vh||1000}});p.on('pageerror',e=>console.log('ERR',e.message));
await p.goto('file://'+path.resolve(__dirname,'..','explorer','LnPc2_Explorer_offline.html'));await p.waitForTimeout(2500);
for(const s of J.jobs){
  const r=await p.evaluate(async s=>{const L=window.LNPC;L.select(s.id);L.clear();L.ui(Object.assign({style:'ball',h:true,geo:false,spin:false,rings:false,spinsph:false,half:false,sym:false,num:false},s.ui||{}));
    if(s.view)L.view(s.view);if(s.focus)L.focus(s.focus,s.fz);(s.rot||[]).forEach(([d,a])=>L.rotate(d,a));if(s.zoom)L.zoom(s.zoom);
    const out={desc:L.desc(),m:[]};if(s.geomOnly){const g=L.geom();return {xyz:g.xyz,el:g.el,m:g.m,url:null}}if(s.lnN)out.lnN=L.lnN();if(s.ln2){const d=L.desc();out.m.push(['dist',[d.m,d.N8[0]],L.meas('dist',[d.m,d.N8[0]])]);out.m.push(['dist',[d.m,d.N8[4]],L.meas('dist',[d.m,d.N8[4]])])}(s.meas||[]).forEach(([t,ix])=>out.m.push([t,ix,L.meas(t,ix)]));
    await new Promise(r=>setTimeout(r,400));const g=L.geom();const add=(a,b)=>a.map((x,i)=>x+b[i]),mul=(a,k)=>a.map(x=>x*k);
    const pts={metal:g.xyz[g.m],Ntop:g.N8.slice(0,4).map(i=>g.xyz[i]),Nbot:g.N8.slice(4).map(i=>g.xyz[i]),Tc:g.T.c,Bc:g.B.c,Tca:add(g.T.c,mul(g.T.v,g.h[0])),Bca:add(g.B.c,mul(g.B.v,g.h[1])),ax:g.ax};
    const flat=[pts.metal,...pts.Ntop,...pts.Nbot,pts.Tc,pts.Bc,pts.Tca,pts.Bca];const rc=document.getElementById('mol').getBoundingClientRect();const pr=L.proj(flat).map(q=>[q[0]-rc.left,q[1]-rc.top]);
    const im=await L.png(s.size||'c2_600',s.ov!==false);out.k=im.k;out.W=im.W;out.H=im.H;out.scr=pr.map(q=>[q[0]*im.k,q[1]*im.k]);out.geo={th:g.th,dN4:g.dN4,dCa:g.dCa};out.url=im.url;return out},s);
  if(r.url)fs.writeFileSync(s.out,Buffer.from(r.url.split(',')[1],'base64'));delete r.url;fs.writeFileSync(s.out.replace('.png','.json'),JSON.stringify(r));console.log(s.out,r.W||'',r.H||'')}
await p.close()}
await b.close()})();
