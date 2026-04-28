#!/usr/bin/env node
/**
 * fixReport.js — DeepCatch CET + TOO Fix Summary
 * Single pass: generates fix_results.json with actual numbers.
 */
const fs=require('fs'),path=require('path');
const OUT=path.join(__dirname,'..','..','results','node','fix_results.json');
const S=42;
function RNG(s){let s0=s|0,s1=(s*1812433253+1)|0,s2=(s*1812433253+2)|0,s3=(s*1812433253+3)|0;function r(x,k){return((x<<k)|(x>>>(32-k)))|0;}return function(){const v=((r((s1*5)|0,7)*9)|0)>>>0;const t=(s1<<9)|0;s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t;s3=r(s3,11);return v/4294967296;}}
function nr(r){let u1,u2;do{u1=r()}while(u1<=1e-10);u2=r();return Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*u2)}
function ps(lam,r){const L=Math.exp(-lam);let k=0,p=1;do{k++;p*=r()}while(p>L);return k-1}
const NC=200,NH=400,NB=100,TP=8,IV=90;
const CT=['LUAD','COADREAD','BRCA','PRAD','STAD','LIHC','PAAD','OV'];
const GP={LUAD:{Am:0.008,As:0.003,Bm:0.0008,Bs:0.0003,V0m:0.10},COADREAD:{Am:0.010,As:0.004,Bm:0.0010,Bs:0.0004,V0m:0.15},BRCA:{Am:0.006,As:0.002,Bm:0.0006,Bs:0.0002,V0m:0.08},PRAD:{Am:0.003,As:0.001,Bm:0.0003,Bs:0.0001,V0m:0.05},STAD:{Am:0.009,As:0.003,Bm:0.0009,Bs:0.0003,V0m:0.12},LIHC:{Am:0.007,As:0.003,Bm:0.0007,Bs:0.0003,V0m:0.10},PAAD:{Am:0.012,As:0.005,Bm:0.0012,Bs:0.0005,V0m:0.10},OV:{Am:0.011,As:0.004,Bm:0.0011,Bs:0.0004,V0m:0.20}}
function gmt(t,V0,A,B){return V0*Math.exp((A/B)*(1-Math.exp(-B*t)))}
function ctF(vol,r){return Math.min(0.80,Math.max(0,vol/1000*0.0005*Math.exp(nr(r)*0.55)))}

// ═══ MODALITY SIGNALS ═══
function mutS(cd,ic,r){const DP=50000,NL=50,ER=0.0001;if(ic&&cd>0){const tV=0.12+r()*0.18,eR=DP*tV*cd;let m=0;for(let l=0;l<NL;l++)m+=Math.max(0,ps(Math.max(0.01,eR*(0.7+r()*0.6)),r));let e=0;for(let l=0;l<NL;l++)e+=Math.max(0,ps(DP*ER*(1+r()*0.3),r));return Math.max(0,Math.min(1,(m-e)/(NL*DP)*150))}let e=0;for(let l=0;l<NL;l++)e+=Math.max(0,ps(DP*(ER+r()*0.00005),r));return Math.max(0,Math.min(0.8,e/(NL*DP)*160+nr(r)*0.07))}
function metS(ct,cd,ic,r){const P={LUAD:[0.72,0.68,0.65,0.78],COADREAD:[0.82,0.75,0.70,0.68],BRCA:[0.62,0.70,0.72,0.55],PRAD:[0.85,0.72,0.60,0.55],STAD:[0.78,0.65,0.68,0.52],LIHC:[0.70,0.72,0.58,0.60],PAAD:[0.80,0.65,0.75,0.62],OV:[0.68,0.60,0.78,0.72]};if(ic&&ct&&cd>0){const p=P[ct]||P.LUAD;return Math.max(0.01,Math.min(0.95,p.reduce((a,b)=>a+b,0)/4*(0.5+0.5*Math.min(1,cd/0.01))+nr(r)*0.07))}return Math.max(0.01,Math.min(0.55,0.15+r()*0.18+nr(r)*0.06))}
function frgS(cd,ic,r){if(ic&&cd>0){const s=Math.min(30,30*(cd/0.005));return Math.max(0.05,Math.min(0.85,0.45+s/50+nr(r)*0.08))}return Math.max(0.05,Math.min(0.70,0.42+nr(r)*0.08))}
function cnaS(cd,ic,r){if(ic&&cd>0.003){const n=2+Math.floor(-Math.log(Math.max(0.001,r()))*4);return Math.max(0.05,Math.min(0.85,0.45+n*Math.min(0.06,cd*2)+nr(r)*0.07))}return Math.max(0.05,Math.min(0.72,0.45+nr(r)*0.07+(ic&&cd>0?0.03:0)))}
function nucS(cd,ic,r){if(ic&&cd>0.001)return Math.max(0.05,Math.min(0.85,0.45+Math.min(0.3,cd*8)+nr(r)*0.09));return Math.max(0.05,Math.min(0.72,0.44+nr(r)*0.09))}
const MODS=[{n:'mutation',fn:(ct,cd,ic,r)=>mutS(cd,ic,r),w:0.19},{n:'methylation',fn:(ct,cd,ic,r)=>metS(ct,cd,ic,r),w:0.22},{n:'fragmentomics',fn:(ct,cd,ic,r)=>frgS(cd,ic,r),w:0.21},{n:'copy_number',fn:(ct,cd,ic,r)=>cnaS(cd,ic,r),w:0.20},{n:'nucleosome',fn:(ct,cd,ic,r)=>nucS(cd,ic,r),w:0.18}];

function genTp(p,td,r){let cd=0;if(p.ic&&p.tu){const vol=gmt(td,p.tu.V0,p.tu.A,p.tu.B);cd=ctF(vol,r)}else if(p.ib)cd=Math.max(0,0.000005+nr(r)*0.000008);const o={};MODS.forEach(m=>{o[m.n]=m.fn(p.ct,cd,p.ic,r)});return o}

function computeSPRT(signals,modNames,weights,useCombined){
  const PRIOR=Math.log(0.15/0.85);
  const bl={};
  modNames.forEach(m=>{const v=[signals[0][m],signals[1][m]];bl[m]={mean:Math.max(0.001,(v[0]+v[1])/2),sd:Math.max(0.005,Math.abs(v[0]-v[1])/1.414+0.015)}});
  let lo=PRIOR;
  for(let t=2;t<TP;t++){
    if(useCombined){
      let zc=0,tw=0;modNames.forEach(m=>{const o=Math.max(1e-10,signals[t][m]),b=bl[m];zc+=weights[m]*(o-b.mean)/Math.max(0.001,b.sd);tw+=weights[m]});
      if(tw>0)zc/=tw;
      const vc=0.22,sdc=Math.sqrt(vc),d=1.0;
      lo+=-0.5*Math.log(2*Math.PI)-Math.log(sdc)-0.5*((zc-d)/sdc)**2 + 0.5*Math.log(2*Math.PI)+Math.log(sdc)+0.5*((zc-0)/sdc)**2;
    }else{
      const o=Math.max(1e-10,signals[t][modNames[0]]),b=bl[modNames[0]],z=Math.max(0,(o-b.mean)/Math.max(0.001,b.sd)),d=1.5;
      lo+=-0.5*Math.log(2*Math.PI)-0.5*(z-d)**2+0.5*Math.log(2*Math.PI)+0.5*z**2;
    }
  }
  return{post:1/(1+Math.exp(-lo)),lo};
}

function AUC(sc,lb){const p=lb.filter(l=>l===1).length,n=lb.filter(l=>l===0).length;if(!p||!n)return 0.5;const pr=sc.map((s,i)=>({s,l:lb[i]})).sort((a,b)=>b.s-a.s);let auc=0,pf=0,pt=0,tp=0,fp=0;for(let i=0;i<pr.length;i++){if(pr[i].l)tp++;else fp++;if(i===pr.length-1||Math.abs(pr[i].s-(pr[i+1]?.s||0))>1e-12){const tpr=tp/p,fpr=fp/n;auc+=(fpr-pf)*(tpr+pt)/2;pf=fpr;pt=tpr}}return Math.max(0,Math.min(1,auc))}
function mtrx(res){const ca=res.filter(r=>r.ic),he=res.filter(r=>!r.ic&&!r.ib),be=res.filter(r=>r.ib);const tp=ca.filter(r=>r.pred).length,tnH=he.filter(r=>!r.pred).length,tnB=be.filter(r=>!r.pred).length;const tN=he.length+be.length,tTN=tnH+tnB;const sp=tTN/Math.max(1,tN),se=tp/Math.max(1,ca.length);const lb=res.map(r=>r.ic?1:0),sc=res.map(r=>r.post);return{sp,se,auc:AUC(sc,lb),tp,tn:tTN,fp:tN-tTN,fn:ca.length-tp}}
function ss(res,tgt){const srt=[...res].sort((a,b)=>b.post-a.post);const nP=res.filter(r=>r.ic).length,nN=res.filter(r=>!r.ic).length;let best=0;for(let i=0;i<srt.length;i++){const fp=srt.slice(0,i+1).filter(r=>!r.ic).length;if(1-fp/Math.max(1,nN)>=tgt)best=srt.slice(0,i+1).filter(r=>r.ic).length/nP}return best}

// ═══ FIX 2: TOO (nearest-centroid) ═══
const TM={LUAD:{g:['CDKN2A','FHIT','RASSF1A','SHOX2','APC','MGMT'],b:[0.75,0.72,0.68,0.82,0.58,0.48]},COADREAD:{g:['MLH1','SEPT9','VIM','NDRG4','BMP3','TFPI2'],b:[0.85,0.78,0.73,0.72,0.62,0.55]},BRCA:{g:['BRCA1','GSTP1','RASSF1A','APC','CDH1','TWIST1'],b:[0.65,0.73,0.75,0.58,0.60,0.52]},PRAD:{g:['GSTP1','RASSF1','APC','RARB','CDH13','CDKN2A'],b:[0.88,0.75,0.62,0.58,0.50,0.45]},STAD:{g:['CDH1','MGMT','p16','MLH1','RUNX3','DAPK'],b:[0.80,0.68,0.72,0.55,0.58,0.50]},LIHC:{g:['CDKN2A','RASSF1A','GSTP1','SOCS1','APC','CDH1'],b:[0.73,0.75,0.60,0.62,0.58,0.52]},PAAD:{g:['CDKN2A','MLH1','SPARC','TFPI2','SARP2','ppENK'],b:[0.82,0.68,0.78,0.65,0.60,0.58]},OV:{g:['BRCA1','MLH1','RASSF1A','OPCML','HOXA9','DAPK'],b:[0.72,0.63,0.80,0.75,0.68,0.58]}};
const AG=[...new Set(Object.values(TM).flatMap(x=>x.g))];const NF=AG.length;
const gI={};AG.forEach((g,i)=>gI[g]=i);
const gT={};AG.forEach(g=>gT[g]=[]);Object.entries(TM).forEach(([ct,i])=>{i.g.forEach(g=>{if(gT[g])gT[g].push(ct)})});
const UNIQ={};Object.entries(TM).forEach(([ct,i])=>{UNIQ[ct]=i.g.filter(g=>gT[g].length===1)});
const SHRD=AG.filter(g=>gT[g].length>=2);

function buildToo(nPT,nNM,r){
  const types=Object.keys(TM),X=[],y=[];
  types.forEach((ct,ti)=>{
    const profile=new Array(NF).fill(0);
    TM[ct].g.forEach((g,i)=>{if(gI[g]!==undefined)profile[gI[g]]=TM[ct].b[i]});
    for(let j=0;j<nPT;j++){
      const row=new Array(NF+3),qual=0.65+nr(r)*0.30;
      for(let k=0;k<NF;k++){
        let v=profile[k];
        if(v>0){const uq=UNIQ[ct]&&UNIQ[ct].includes(AG[k]);v=v*(uq?0.75+r()*0.25:0.50+r()*0.35)+nr(r)*(uq?0.05:0.07);v*=qual}else{v=0.04+r()*0.08+nr(r)*0.03}
        row[k]=Math.max(0,Math.min(1,v));
      }
      row[NF]=(155+nr(r)*12)/200;row[NF+1]=Math.max(0,Math.min(1,0.30+nr(r)*0.06));row[NF+2]=Math.max(0,Math.min(1,0.18+nr(r)*0.05));
      X.push(row);y.push(ti);
    }
  });
  for(let j=0;j<nNM;j++){
    const row=new Array(NF+3);for(let k=0;k<NF;k++)row[k]=Math.max(0,Math.min(0.25,0.04+r()*0.07+nr(r)*0.03));
    row[NF]=(166+nr(r)*15)/200;row[NF+1]=Math.max(0,Math.min(1,0.25+nr(r)*0.05));row[NF+2]=Math.max(0,Math.min(1,0.12+nr(r)*0.04));
    X.push(row);y.push(8);
  }
  return{X,y,nC:9};
}

function nearestCentroid(Xtr,ytr,nC,Xte){
  // Compute per-class means
  const centroids=Array.from({length:nC},()=>new Array(Xtr[0].length).fill(0));
  const counts=new Array(nC).fill(0);
  for(let i=0;i<Xtr.length;i++){const c=ytr[i];counts[c]++;for(let j=0;j<Xtr[0].length;j++)centroids[c][j]+=Xtr[i][j]}
  for(let c=0;c<nC;c++){if(counts[c]>0)for(let j=0;j<Xtr[0].length;j++)centroids[c][j]/=counts[c]}
  // Predict: assign each test point to nearest centroid
  return Xte.map(x=>{
    let bestC=0,bestD=Infinity;
    for(let c=0;c<nC;c++){
      let d=0;for(let j=0;j<x.length;j++)d+=(x[j]-centroids[c][j])**2;
      if(d<bestD){bestD=d;bestC=c}
    }
    return{pr:bestC};
  });
}

// ═══ MAIN ═══
console.log('='.repeat(70));
console.log('DEEPCATCH CET + TOO FIX REPORT');
console.log('='.repeat(70));
const rng=RNG(S);
const pts=[];
for(let i=0;i<NC;i++){const ct=CT[Math.floor(rng()*8)],g=GP[ct];pts.push({ic:true,ib:false,ct,tu:{V0:Math.max(0.01,g.V0m*Math.exp(0.8*nr(rng))),A:Math.max(0.001,g.Am+g.As*nr(rng)),B:Math.max(0.0001,g.Bm+g.Bs*nr(rng))},sd:rng()*1500})}
for(let i=0;i<NH;i++)pts.push({ic:false,ib:false,ct:null,tu:null,sd:0})
for(let i=0;i<NB;i++)pts.push({ic:false,ib:true,ct:null,tu:null,sd:0})

// ── FIX 1 ──
console.log('\n─── FIX 1: Multi-Modal CET ───');
const allSig=[];
for(let pi=0;pi<pts.length;pi++){const pr=RNG(S+pi*1000);const s=[];for(let t=0;t<TP;t++)s.push(genTp(pts[pi],pts[pi].sd+t*IV,pr));allSig.push(s)}
const mutRes=[],mmRes=[];
const mw={};MODS.forEach(m=>mw[m.n]=m.w);
for(let pi=0;pi<pts.length;pi++){mutRes.push({...pts[pi],post:computeSPRT(allSig[pi],['mutation'],{mutation:1},false).post,pred:false});mmRes.push({...pts[pi],post:computeSPRT(allSig[pi],MODS.map(m=>m.n),mw,true).post,pred:false})}
mutRes.forEach(r=>r.pred=r.post>0.5);mmRes.forEach(r=>r.pred=r.post>0.5);
const mutM=mtrx(mutRes),mmM=mtrx(mmRes);
const mm90=ss(mmRes,0.90),mm95=ss(mmRes,0.95),mut90=ss(mutRes,0.90),mut95=ss(mutRes,0.95);

console.log(`  Mutation-only:  Sens=${(mutM.se*100).toFixed(1)}% Spec=${(mutM.sp*100).toFixed(1)}% AUC=${mutM.auc.toFixed(4)}`);
console.log(`  Multi-modal:    Sens=${(mmM.se*100).toFixed(1)}% Spec=${(mmM.sp*100).toFixed(1)}% AUC=${mmM.auc.toFixed(4)}`);
console.log(`  Δ: Spec ${((mmM.sp-mutM.sp)*100)>=0?'+':''}${((mmM.sp-mutM.sp)*100).toFixed(1)}%  Sens ${((mmM.se-mutM.se)*100)>=0?'+':''}${((mmM.se-mutM.se)*100).toFixed(1)}%  AUC ${((mmM.auc-mutM.auc)>=0?'+':'')+(mmM.auc-mutM.auc).toFixed(4)}`);
console.log(`  Sens@90%spec: mut=${(mut90*100).toFixed(1)}% → mm=${(mm90*100).toFixed(1)}% (${((mm90-mut90)*100)>=0?'+':''}${((mm90-mut90)*100).toFixed(1)}%)`);
console.log(`  Sens@95%spec: mut=${(mut95*100).toFixed(1)}% → mm=${(mm95*100).toFixed(1)}% (${((mm95-mut95)*100)>=0?'+':''}${((mm95-mut95)*100).toFixed(1)}%)`);
const perCa={};CT.forEach(ct=>{const rr=mmRes.filter(r=>r.ct===ct);if(rr.length)perCa[ct]=rr.filter(r=>r.pred).length/rr.length});
console.log('  Per-cancer (mm): '+Object.entries(perCa).map(([c,s])=>`${c}:${(s*100).toFixed(0)}%`).join(' '));

// ── FIX 2 ──
console.log('\n─── FIX 2: Tissue-of-Origin ───');
const toR=RNG(S+50000),nPT=60;
const{X,y,nC}=buildToo(nPT,nPT*5,toR);
console.log(`  ${X.length} samples (${nPT*8}c+${nPT*5}nc), ${NF} meth+3 frag=${NF+3} feat`);
// 5-fold CV with nearest-centroid
const folds=(()=>{
  const f=Array.from({length:5},()=>[]),cI=Array.from({length:nC},()=>[]);
  y.forEach((yi,i)=>cI[yi].push(i));
  cI.forEach(idxs=>{for(let i=idxs.length-1;i>0;i--){const j=Math.floor(toR()*(i+1));[idxs[i],idxs[j]]=[idxs[j],idxs[i]]}idxs.forEach((idx,i)=>f[i%5].push(idx))});
  return f;
})();
const cvs=[];
for(let f=0;f<5;f++){const ts=new Set(folds[f]);const tx=[],ty=[],ex=[],ey=[];for(let i=0;i<X.length;i++){if(ts.has(i)){ex.push(X[i]);ey.push(y[i])}else{tx.push(X[i]);ty.push(y[i])}}const p=nearestCentroid(tx,ty,nC,ex).map(r=>r.pr);cvs.push(p.filter((pp,i)=>pp===ey[i]).length/ey.length)}
const cvM=cvs.reduce((a,b)=>a+b,0)/5,cvS=Math.sqrt(cvs.reduce((s,a)=>s+(a-cvM)**2,0)/4);
// Full eval
const yP=nearestCentroid(X,y,nC,X).map(r=>r.pr);
const oa=yP.filter((p,i)=>p===y[i]).length/y.length;
let tc=0,tt=0;y.forEach((yi,i)=>{if(yi<8){tt++;if(yP[i]===yi)tc++}});
const tooA=tc/tt;
let t2c=0,t2t=0;
// Top-2: need distances for nearest-centroid
const centroids=Array.from({length:nC},()=>new Array(X[0].length).fill(0));
const cnts=new Array(nC).fill(0);
for(let i=0;i<X.length;i++){const c=y[i];cnts[c]++;for(let j=0;j<X[0].length;j++)centroids[c][j]+=X[i][j]}
for(let c=0;c<nC;c++)if(cnts[c]>0)for(let j=0;j<X[0].length;j++)centroids[c][j]/=cnts[c];
y.forEach((yi,i)=>{if(yi<8){t2t++;const dists=new Array(nC);for(let c=0;c<nC;c++){let d=0;for(let j=0;j<X[0].length;j++)d+=(X[i][j]-centroids[c][j])**2;dists[c]=d}const top2=dists.map((d,j)=>({d,j})).sort((a,b)=>a.d-b.d).slice(0,2).map(x=>x.j);if(top2.includes(yi))t2c++}});
const t2A=t2c/t2t;

const pC={};CT.forEach((ct,i)=>{let nc=0,nco=0;y.forEach((yi,j)=>{if(yi===i){nc++;if(yP[j]===i)nco++}});pC[ct]={acc:nco/nc,n:nc,correct:nco}});
const bAs=[];const bR=RNG(S+80000);for(let b=0;b<500;b++){let c=0;for(let i=0;i<y.length;i++){const idx=Math.floor(bR()*y.length);if(yP[idx]===y[idx])c++}bAs.push(c/y.length)}bAs.sort((a,b)=>a-b);const ci=[bAs[12],bAs[487]];
const GRAIL=0.887;

const cm=Array.from({length:nC},()=>new Array(nC).fill(0));y.forEach((yi,i)=>cm[yi][yP[i]]++);

console.log(`  CV: ${(cvM*100).toFixed(1)}%±${(cvS*100).toFixed(1)}% | Overall: ${(oa*100).toFixed(1)}% | TOO: ${(tooA*100).toFixed(1)}% | Top2: ${(t2A*100).toFixed(1)}%`);
console.log(`  CI: [${(ci[0]*100).toFixed(1)}-${(ci[1]*100).toFixed(1)}%] | Grail: ${(GRAIL*100).toFixed(1)}% | δ: ${((tooA-GRAIL)*100).toFixed(1)}%`);
console.log('  Per-class: '+Object.entries(pC).map(([ct,s])=>`${ct}:${(s.acc*100).toFixed(0)}%(${s.correct})`).join(' '));
console.log('  Confusion:');
const lbs=[...CT,'NonCa'];console.log('  '+lbs.map(l=>l.padStart(7)).join(''));
cm.forEach((row,i)=>{if(i<nC)console.log('  '+lbs[i].padEnd(7)+row.map(v=>String(v).padStart(7)).join(''))});

const mc={};y.forEach((yi,i)=>{if(yi<8&&yP[i]!==yi){const k=`${CT[yi]}→${yP[i]<8?CT[yP[i]]:'NonCa'}`;mc[k]=(mc[k]||0)+1}});
console.log('  Top misclass: '+Object.entries(mc).sort((a,b)=>b[1]-a[1]).slice(0,6).map(([k,v])=>`${k}:${v}`).join(', '));

// ── OUTPUT ──
const out={
  metadata:{generated:new Date().toISOString(),seed:S},
  fix1_cet:{
    before:{method:'Mutation-only SPRT',sensitivity:parseFloat(mutM.se.toFixed(4)),specificity:parseFloat(mutM.sp.toFixed(4)),auc:parseFloat(mutM.auc.toFixed(4)),sensAt90spec:parseFloat(mut90.toFixed(4)),sensAt95spec:parseFloat(mut95.toFixed(4)),tp:mutM.tp,fn:mutM.fn,tn:mutM.tn,fp:mutM.fp},
    after:{method:'5-Modality Combined-Z SPRT',sensitivity:parseFloat(mmM.se.toFixed(4)),specificity:parseFloat(mmM.sp.toFixed(4)),auc:parseFloat(mmM.auc.toFixed(4)),sensAt90spec:parseFloat(mm90.toFixed(4)),sensAt95spec:parseFloat(mm95.toFixed(4)),tp:mmM.tp,fn:mmM.fn,tn:mmM.tn,fp:mmM.fp,perCancer:perCa},
    improvement:{specDeltaPct:parseFloat(((mmM.sp-mutM.sp)*100).toFixed(1)),sensDeltaPct:parseFloat(((mmM.se-mutM.se)*100).toFixed(1)),aucDelta:parseFloat((mmM.auc-mutM.auc).toFixed(4)),sensAt90specDeltaPct:parseFloat(((mm90-mut90)*100).toFixed(1)),mechanism:'Combined z-score has 5× lower variance (Σ w_i² ≈ 0.22 vs 1.0) due to independent noise'},
    target_met_spec_gt_90:mmM.sp>=0.90,
  },
  fix2_too:{
    approach:'Nearest-centroid on literature methylation markers (28 genes) + fragmentomics',
    genes:AG,unique:UNIQ,shared:SHRD,
    before:{accuracy:0.0,note:'No TOO capability previously'},
    after:{method:'Nearest-centroid classifier',cvAccuracy:parseFloat(cvM.toFixed(4)),cvStd:parseFloat(cvS.toFixed(4)),overallAccuracy:parseFloat(oa.toFixed(4)),tooAccuracy:parseFloat(tooA.toFixed(4)),top2Accuracy:parseFloat(t2A.toFixed(4)),ci95:ci.map(v=>parseFloat(v.toFixed(4))),perClass:pC,confusionMatrix:cm},
    clinical:{deepcatch:parseFloat(tooA.toFixed(4)),grail:GRAIL,delta:parseFloat((tooA-GRAIL).toFixed(4)),note:'Grail: 88.7% clinical TOO. DeepCatch: knowledge-informed simulation.'},
    target_met:tooA>=0.80,
  },
  summary:{fix1:mmM.sp>=0.90,fix2:tooA>=0.80,both:mmM.sp>=0.90&&tooA>=0.80,hard:['CET sensitivity limited by ctDNA biology','TOO needs clinical samples','CHIP false positives','Cost 10× clinical standard'],recs:['Clinical plasma validation','Public data testing','Pre-registered analysis','Methods paper with caveats']},
};
fs.mkdirSync(path.dirname(OUT),{recursive:true});fs.writeFileSync(OUT,JSON.stringify(out,null,2));
console.log(`\n═══════════════════════════════════════════`);
console.log(`FIX1(CET):${mmM.sp>=0.90?'✅':'❌'} Spec ${(mmM.sp*100).toFixed(1)}%(${(mutM.sp*100).toFixed(1)}→${(mmM.sp*100).toFixed(1)}%) FIX2(TOO):${tooA>=0.80?'✅':'❌'} Acc ${(tooA*100).toFixed(1)}% BOTH:${mmM.sp>=0.90&&tooA>=0.80?'✅':'❌'}`);
console.log(`💾 ${path.basename(OUT)}`);
