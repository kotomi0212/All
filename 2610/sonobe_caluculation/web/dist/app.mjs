import {readRecording,analyze,AXES,demean,resultCSV} from './analysis.mjs';
const $=id=>document.getElementById(id);
let recording=null,result=null,axis='AP',loadId=0;
const format=(n,d=3)=>n==null?'—':n.toFixed(d);
function status(text,kind=''){$('analysisstatus').textContent=text;$('analysisstatus').className='status '+kind;}
function message(text){$('message').textContent=text;$('message').hidden=!text;}
function resetResults(){result=null;$('results').hidden=true;$('empty').hidden=false;}
function controls(enabled){for(const id of ['start','end','recalc','full'])$(id).disabled=!enabled;}
async function loadFile(file) {
  if(!file)return;
  const id=++loadId;recording=null;resetResults();controls(false);message('');
  $('filename').textContent=file.name;$('filemeta').textContent='ファイルを読み込んでいます…';status('解析中');
  try {
    if(file.size>50*1024*1024)throw new Error('ファイルが50 MBを超えています。解析する記録を分けてください。');
    const buffer=await file.arrayBuffer();if(id!==loadId)return;
    let text;try{text=new TextDecoder('utf-8',{fatal:true}).decode(buffer);}
    catch{text=new TextDecoder('shift_jis',{fatal:true}).decode(buffer);}
    recording=readRecording(text,file.name);
    const r=recording;
    $('filemeta').textContent=`${format(r.fs,0)} Hz · ${r.time.length.toLocaleString()}点 · ${format(r.max-r.min,2)}秒`;
    const standard=r.min<=5&&r.max>=35-1e-7;
    $('start').value=standard?5:r.min;$('end').value=standard?35:Number(r.max.toFixed(8));
    $('intervalnote').textContent=standard?'初期設定：5〜35秒（論文の解析区間）':'5〜35秒が範囲外のため、全区間を設定しました。';
    controls(true);calculate();
  } catch(e) {recording=null;resetResults();controls(false);message(e.message);status('読み込みエラー');$('filemeta').textContent='ファイル形式を確認してください。';}
}

function calculate() {
  if(!recording)return;
  try {
    if(!$('start').value.trim()||!$('end').value.trim())throw new Error('開始時刻と終了時刻を入力してください。');
    result=analyze(recording,Number($('start').value),Number($('end').value));
    message('');$('empty').hidden=true;$('results').hidden=false;
    status(`${format(result.start_s,2)}–${format(result.end_s,2)} s · ${result.samples.toLocaleString()}点`,'ready');
    renderCards();renderWarnings();renderMetrics();renderChart();
  } catch(e){resetResults();status('区間を確認');message(e.message);}
}

function renderCards() {
  $('cards').replaceChildren();
  for(const key of ['AP','ML']) {
    const m=result.axes[key],card=document.createElement('article');card.className='mode-card';
    // All interpolated values below come from validated numeric calculations or fixed labels.
    card.innerHTML=`<div class="card-head"><h3>${key==='AP'?'前後方向':'左右方向'}</h3><span>${key}</span></div><div class="percentages"><div class="percentage">足関節<strong>${format(m.ankle_percent,1)}<small>%</small></strong></div><div class="percentage">股関節<strong>${format(m.hip_percent,1)}<small>%</small></strong></div></div><div class="share-bar" role="img" aria-label="足関節 ${format(m.ankle_percent,1)}%、股関節 ${format(m.hip_percent,1)}%"><span style="width:${m.ankle_percent??0}%"></span><span style="width:${m.hip_percent??0}%"></span></div>`;
    $('cards').append(card);
  }
}
function renderWarnings() {
  const box=$('warnings');box.replaceChildren();box.hidden=!result.warnings.length;
  if(result.warnings.length){const h=document.createElement('strong');h.textContent='データの整合性を確認してください';box.append(h);for(const w of result.warnings){const p=document.createElement('p');p.textContent=w;box.append(p);}}
}
function renderMetrics() {
  const entries=[['足関節 RMS (mm/s²)',m=>format(m.ankle_rms)],['股関節 RMS (mm/s²)',m=>format(m.hip_rms)],['RXI（股 / 足）',m=>format(m.rxi)],['CXI（相関係数）',m=>format(m.cxi)],['モード和と頭部のRMSE (mm/s²)',m=>format(m.reconstruction_rmse)],['整合性（確認目安：頭部RMSの1%）',m=>m.consistent===null?'未確認':m.consistent?'目安の範囲内':'要確認']];
  $('metrics').replaceChildren();
  for(const [label,fn] of entries){const tr=document.createElement('tr');for(const text of [label,fn(result.axes.AP),fn(result.axes.ML)]){const td=document.createElement('td');td.textContent=text;if(text==='要確認')td.className='check';tr.append(td);}$('metrics').append(tr);}
}

function chartSVG() {
  const r=recording, ids=result.ids,columns=AXES[axis].filter(c=>r.signals[c]);
  const series=columns.map(c=>demean(ids.map(i=>r.signals[c][i])));
  let max=0;for(const a of series)for(const v of a)max=Math.max(max,Math.abs(v));
  max=max||1;
  const exp=10**Math.floor(Math.log10(max)),top=Math.ceil(max/exp)*exp;
  const W=800,H=300,L=65,R=22,T=28,B=46,pw=W-L-R,ph=H-T-B;
  const x=t=>L+(t-result.start_s)/(result.end_s-result.start_s)*pw;
  const y=v=>T+(top-v)/(2*top)*ph;
  let body=`<rect width="${W}" height="${H}" fill="white"/><text x="${L}" y="14" fill="#7a887e" font-size="10">${axis} · mm/s²</text>`;
  for(let i=0;i<=4;i++) {const v=top-i*top/2;body+=`<line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="#e8ece7"/><text x="${L-10}" y="${y(v)+4}" text-anchor="end" font-size="11" fill="#7a887e">${Number(v.toPrecision(4))}</text>`;}
  for(let i=0;i<=6;i++){const t=result.start_s+i*(result.end_s-result.start_s)/6;body+=`<line x1="${x(t)}" x2="${x(t)}" y1="${T}" y2="${H-B}" stroke="#f0f3ed"/><text x="${x(t)}" y="${H-B+22}" text-anchor="middle" fill="#7a887e" font-size="11">${Number(t.toFixed(2))}</text>`;}
  const colors=['#527d6c','#c38667','#7f8f9c'];
  for(let j=0;j<series.length;j++) {
    // Min/max decimation retains peaks. All samples still enter the statistics.
    const a=series[j], points=[],bucket=Math.max(1,Math.ceil(a.length/900));
    for(let i=0;i<a.length;i+=bucket){let lo=i,hi=i;for(let k=i+1;k<Math.min(i+bucket,a.length);k++){if(a[k]<a[lo])lo=k;if(a[k]>a[hi])hi=k;}for(const k of [...new Set([i,lo,hi,Math.min(i+bucket-1,a.length-1)])].sort((u,v)=>u-v))points.push(`${x(r.time[ids[k]]).toFixed(2)},${y(a[k]).toFixed(2)}`);}
    body+=`<polyline points="${points.join(' ')}" fill="none" stroke="${colors[j]}" stroke-width="1.2" opacity=".85"/>`;
  }
  body+=`<text x="${W/2}" y="${H-5}" text-anchor="middle" fill="#7a887e" font-size="11">時刻 (s)</text>`;
  const warning=result.warnings.filter(w=>w.startsWith(axis+':')).join(' ');
  const escape=s=>s.replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('>','&gt;');
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-labelledby="charttitle chartdesc" font-family="sans-serif"><title id="charttitle">${axis} 足関節・股関節・頭部の加速度波形</title><desc id="chartdesc">足関節は緑、股関節は茶、頭部は灰色。${escape(warning)} 選択区間の平均を除去。表示はピークを保持して間引き。</desc>${body}</svg>`;
}
function renderChart(){if(!result)return;$('chart').innerHTML=chartSVG();$('chartcaption').textContent=`${result.samples.toLocaleString()}点を計算に使用。表示のみ、ピークを保持して間引いています。`;}
function download(content,type,suffix) {
  const blob=new Blob([content],{type}),url=URL.createObjectURL(blob),a=document.createElement('a');
  a.href=url;a.download=`balance_${axis}_${result.start_s}-${result.end_s}.${suffix}`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
$('file').addEventListener('change',e=>loadFile(e.target.files[0]));
for(const name of ['dragenter','dragover'])$('dropzone').addEventListener(name,e=>{e.preventDefault();$('dropzone').classList.add('dragover');});
for(const name of ['dragleave','drop'])$('dropzone').addEventListener(name,e=>{e.preventDefault();$('dropzone').classList.remove('dragover');});
$('dropzone').addEventListener('drop',e=>{const files=e.dataTransfer.files;if(files.length>1){message('一度に1つのCSVを選んでください。');return;}loadFile(files[0]);});
document.addEventListener('dragover',e=>e.preventDefault());document.addEventListener('drop',e=>e.preventDefault());
$('rangeform').addEventListener('submit',e=>{e.preventDefault();calculate();});
for(const id of ['start','end'])$(id).addEventListener('input',()=>{if(result)status('区間変更あり · 再解析してください','pending');});
$('full').addEventListener('click',()=>{if(recording){$('start').value=recording.min;$('end').value=Number(recording.max.toFixed(8));calculate();}});
for(const key of ['AP','ML'])$(key.toLowerCase()).addEventListener('click',()=>{axis=key;for(const k of ['AP','ML'])$(k.toLowerCase()).setAttribute('aria-pressed',String(k===axis));renderChart();});
$('csvexport').addEventListener('click',()=>{if(result)download(resultCSV(result),'text/csv;charset=utf-8','csv');});
$('jsonexport').addEventListener('click',()=>{if(result){const{ids,...data}=result;download(JSON.stringify(data,null,2),'application/json','json');}});
$('svgexport').addEventListener('click',()=>{if(result)download(chartSVG(),'image/svg+xml','svg');});
$('clear').addEventListener('click',()=>{loadId++;recording=null;resetResults();controls(false);$('file').value='';$('filename').textContent='ファイル未選択';$('filemeta').textContent='読み込んだデータを消去しました。';$('start').value=5;$('end').value=35;$('intervalnote').textContent='初期設定：5〜35秒（論文の解析区間）';message('');status('データ待ち');});
