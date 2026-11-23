(() => {
 'use strict';
 const TwoLink=globalThis.TwoLink;
 if(!TwoLink||typeof TwoLink.createSimulation!=='function')throw new Error('身体モデルを読み込めませんでした。最新版の two-link-lab.html を開いてください。');
 const $=id=>{const element=document.getElementById(id);if(!element)throw new Error('画面の部品が見つかりません: '+id);return element},P=TwoLink.params,deg=180/Math.PI,dt=1/240;
 let s=Array(8).fill(0),t=0,phase=0,acc=0,prev=0,auto=false,paused=false,stopped=false,drag=null,history=[],sample=0;
 let simulation=TwoLink.createSimulation(),reference=TwoLink.createSimulation(),sRef=Array(8).fill(0);
 const input={floor:0,scene:0,control:true,weight:.5};
 function labels(){const visual=$('mode').value==='visual';$('position-value').value=(+$('position').value/100*(visual?3:8)).toFixed(1)+(visual?'°':' cm');$('frequency-value').value=(+$('frequency').value).toFixed(2)+' Hz';$('amplitude-value').value=(+$('amplitude').value).toFixed(1)+(visual?'°':' cm');$('weight-value').value=(+$('weight').value).toFixed(2);$('hint').textContent=visual?'左右にドラッグ → 視覚背景を傾ける（±3°）':'左右にドラッグ → 床を移動（±8 cm）'}
 function profileLabels(){const profile=TwoLink.profiles[$('profile').value];$('selected-title').textContent=profile.label;$('profile-note').textContent=$('profile').value==='reference'?'左右とも従来の教育用モデルです。高齢者参考条件に切り替えると、文献の群間差による感度分析を行います。':`右側の感度分析：K × ${profile.kScale.toFixed(3)}、追加遅延 ${(profile.extraDelay*1000).toFixed(0)} ms。人体の年齢別同定モデルではありません。`}
 function reset(){s=Array(8).fill(0);sRef=Array(8).fill(0);simulation=TwoLink.createSimulation($('profile').value);reference=TwoLink.createSimulation();t=phase=acc=sample=0;history=[];stopped=false;paused=false;drag=null;$('position').value=0;$('pause').textContent='一時停止';labels();profileLabels()}
 function source(a){auto=a;$('manual').classList.toggle('active',!a);$('auto').classList.toggle('active',a)}
 $('manual').onclick=()=>source(false);$('auto').onclick=()=>source(true);$('reset').onclick=reset;
 $('pause').onclick=()=>{if(stopped)return;paused=!paused;$('pause').textContent=paused?'再開':'一時停止'};
 $('mode').onchange=reset;
 $('profile').onchange=reset;
 ['position','frequency','amplitude','weight'].forEach(id=>$(id).oninput=()=>{if(id==='position')source(false);labels()});
 document.querySelectorAll('[data-f]').forEach(b=>b.onclick=()=>{$('frequency').value=b.dataset.f;source(true);labels()});
 const world=$('world'),referenceWorld=$('world-reference');
 [world,referenceWorld].forEach(canvas=>{
  canvas.onpointerdown=e=>{source(false);drag={x:e.clientX,start:+$('position').value};canvas.setPointerCapture(e.pointerId)};
  canvas.onpointermove=e=>{if(!drag)return;$('position').value=Math.max(-100,Math.min(100,drag.start+(e.clientX-drag.x)/canvas.getBoundingClientRect().width*240));labels()};
  canvas.onpointerup=canvas.onpointercancel=canvas.onlostpointercapture=()=>drag=null;
 });
 document.addEventListener('visibilitychange',()=>{prev=0;acc=0});
 function size(canvas){const r=canvas.getBoundingClientRect(),d=Math.min(globalThis.devicePixelRatio||1,2);if(canvas.width!==Math.round(r.width*d)||canvas.height!==Math.round(r.height*d)){canvas.width=Math.round(r.width*d);canvas.height=Math.round(r.height*d)}const c=canvas.getContext('2d');if(!c)throw new Error('このブラウザではCanvas 2Dを利用できません。');c.setTransform(d,0,0,d,0,0);c.clearRect(0,0,r.width,r.height);return {c,w:r.width,h:r.height}}
 function line(c,points,color,width=2){c.strokeStyle=color;c.lineWidth=width;c.lineCap='round';c.beginPath();points.forEach(([x,y],i)=>i?c.lineTo(x,y):c.moveTo(x,y));c.stroke()}
 function drawBody(canvas,body){
  const {c,w,h}=size(canvas),commonWidth=Math.min(world.getBoundingClientRect().width,referenceWorld.getBoundingClientRect().width),scale=Math.min(h*.57,commonWidth*.9),ground=h*.83,ox=w/2+body[4]*scale;
  c.save();c.translate(w/2,ground);c.rotate(body[6]);for(let x=-w;x<w;x+=48)line(c,[[x,-h],[x,30]],'#26384b',1);c.restore();
  c.setLineDash([4,6]);line(c,[[w/2,ground],[w/2,45]],'#748596',1);c.setLineDash([]);
  const hip=[ox+Math.sin(body[0])*P.l1*scale,ground-Math.cos(body[0])*P.l1*scale],head=[hip[0]+Math.sin(body[1])*P.l2*scale,hip[1]-Math.cos(body[1])*P.l2*scale];
  c.fillStyle='#385365';c.fillRect(0,ground+8,w,25);for(let i=-15;i<16;i++){const x=w/2+i*40+body[4]*scale;line(c,[[x,ground+9],[x-14,ground+31]],'#6c8593',2)}
  c.fillStyle='#d9e4ee';c.fillRect(ox-20,ground-3,43,10);line(c,[[ox,ground],hip],'#59ddc8',15);line(c,[hip,head],'#ffab79',20);
  c.fillStyle='#e8eff5';for(const p of [[ox,ground],hip]){c.beginPath();c.arc(...p,6,0,Math.PI*2);c.fill()}c.beginPath();c.arc(head[0],head[1]-13,13,0,Math.PI*2);c.fill();
  const comX=body[4]+(P.H[0]*body[0]+P.H[1]*body[1])/(P.m1+P.m2),comY=(P.H[0]+P.H[1])/(P.m1+P.m2);
  c.fillStyle='#f2d777';c.beginPath();c.arc(w/2+comX*scale,ground-comY*scale,5,0,Math.PI*2);c.fill();
  c.font='12px system-ui';c.fillStyle='#aebdce';c.fillText('● 重心',12,h-18);c.fillText('背景 '+(body[6]*deg).toFixed(2)+'°',12,h-39);line(c,[[w-78,h-28],[w-78+.1*scale,h-28]],'#aebdce',2);c.fillText('10 cm',w-78,h-38);
 }
 function draw(){
  drawBody(referenceWorld,sRef);drawBody(world,s);
  $('leg').textContent=(s[0]*deg).toFixed(2)+'°';$('hat').textContent=(s[1]*deg).toFixed(2)+'°';$('hip').textContent=((s[1]-s[0])*deg).toFixed(2)+'°';$('floor').textContent=(s[4]*100).toFixed(1)+' cm';
  $('leg-reference').textContent=(sRef[0]*deg).toFixed(2)+'°';$('hat-reference').textContent=(sRef[1]*deg).toFixed(2)+'°';$('hip-reference').textContent=((sRef[1]-sRef[0])*deg).toFixed(2)+'°';$('floor-reference').textContent=(sRef[4]*100).toFixed(1)+' cm';
  $('status').textContent=stopped?'いずれかの条件が20°に達したため停止。リセットしてください。':paused?'一時停止中':(auto?'正弦波入力':'手動入力')+' / '+(input.control?'固定の姿勢制御あり':'制御なし・不安定');
  const showRef=$('overlay').checked&&$('profile').value!=='reference';
  const chart=size($('trace')),cx=chart.c,cw=chart.w,ch=chart.h,mx=Math.max(1,...history.flatMap(v=>v.slice(1,showRef?5:4).map(Math.abs)))*1.15;
  cx.font='11px system-ui';for(const v of [-mx,0,mx]){const y=ch/2-v/mx*(ch/2-18);line(cx,[[42,y],[cw-10,y]],'#334353',1);cx.fillStyle='#9cacc0';cx.fillText(v.toFixed(1)+'°',2,y+4)}
  ['#59ddc8','#ffab79','#c6b5ff'].forEach((col,j)=>line(cx,history.map(v=>[42+(v[0]-Math.max(0,t-12))/12*(cw-52),ch/2-v[j+1]/mx*(ch/2-18)]),col));
  if(showRef){cx.setLineDash([5,4]);line(cx,history.map(v=>[42+(v[0]-Math.max(0,t-12))/12*(cw-52),ch/2-v[4]/mx*(ch/2-18)]),'#f2f5fa');cx.setLineDash([])}
  if(history.length){const rms=j=>Math.sqrt(history.reduce((sum,v)=>sum+v[j]**2,0)/history.length);$('comparison').textContent=`${TwoLink.profiles[$('profile').value].label} / 表示区間の股関節角RMS：選択条件 ${rms(3).toFixed(3)}°・基準 ${rms(4).toFixed(3)}°（過渡応答を含む）`}else $('comparison').textContent='条件を選び、同じ入力で基準との差を観察できます。';
 }
 function frame(now){if(!prev)prev=now;const elapsed=Math.min((now-prev)/1000,.08);prev=now;
  if(!paused&&!stopped){acc+=elapsed;while(acc>=dt){const visual=$('mode').value==='visual';phase+=2*Math.PI*(+$('frequency').value)*dt;const value=auto?Math.sin(phase)*(+$('amplitude').value):+$('position').value/100*(visual?3:8);input.floor=visual?0:value/100;input.scene=visual?value/deg:0;input.control=$('control').checked;input.weight=+$('weight').value;s=simulation.step(input,dt);sRef=reference.step(input,dt);t+=dt;acc-=dt;if(++sample%8===0){history.push([t,s[0]*deg,s[1]*deg,(s[1]-s[0])*deg,(sRef[1]-sRef[0])*deg]);while(history.length&&history[0][0]<t-12)history.shift()}if([...s,...sRef].some(v=>!Number.isFinite(v))||Math.max(...s.slice(0,2).map(Math.abs),...sRef.slice(0,2).map(Math.abs))>20/deg){stopped=true;acc=0;break}}}else acc=0;
  draw();if(globalThis.twoLinkBootComplete)globalThis.twoLinkBootComplete();requestAnimationFrame(frame);
 }
 $('parameters').innerHTML='<table><tr><th>量</th><th>脚（両脚）</th><th>HAT（頭・腕・体幹）</th></tr>'+[['質量 / kg',P.m1,P.m2],['長さ / m',P.l1,P.l2],['近位端から重心 / m',P.h1,P.h2],['近位端まわり慣性 / kg·m²',P.J1,P.J2]].map(r=>'<tr><td>'+r[0]+'</td><td>'+r[1].toFixed(3)+'</td><td>'+r[2].toFixed(3)+'</td></tr>').join('')+'</table>';
 reset();requestAnimationFrame(frame);
})();

