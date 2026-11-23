export const REFERENCE = 'https://doi.org/10.3389/fneur.2025.1500308';
export const AXES = {AP:['xix1','xix2','axhd'],ML:['xiy1','xiy2','ayhd']};
const MODES=['xix1','xix2','xiy1','xiy2'];

export function parseCSV(text) {
  const rows=[]; let row=[], field='', quoted=false;
  for(let i=0;i<text.length;i++) {
    const c=text[i];
    if(quoted) {
      if(c==='"' && text[i+1]==='"'){field+='"';i++;}
      else if(c==='"') quoted=false;
      else field+=c;
    } else if(c==='"' && field==='') quoted=true;
    else if(c===','){row.push(field);field='';}
    else if(c==='\n'||c==='\r') {
      if(c==='\r'&&text[i+1]==='\n')i++;
      row.push(field);rows.push(row);row=[];field='';
    } else field+=c;
  }
  if(quoted)throw new Error('CSVの引用符が閉じられていません。');
  if(field!==''||row.length){row.push(field);rows.push(row);}
  return rows;
}

export function readRecording(text,name='recording.csv') {
  const rows=parseCSV(text.replace(/^\uFEFF/,''));
  const heads=[];
  rows.forEach((r,i)=>{if(MODES.every(c=>r.includes(c)))heads.push(i);});
  if(!heads.length)throw new Error('処理済みモード列（xix1, xix2, xiy1, xiy2）が見つかりません。MedicalForcePlateの解析結果付きCSVを選んでください。生の床反力・IMUのみのCSVには未対応です。');
  const hi=heads.at(-1), header=rows[hi].map(x=>x.trim());
  if(new Set(header).size!==header.length)throw new Error('時系列の列名が重複しています。');
  const names=[...MODES,...['axhd','ayhd'].filter(c=>header.includes(c))];
  const indexes=names.map(c=>header.indexOf(c));
  const signals=Object.fromEntries(names.map(n=>[n,[]])), time=[];
  for(let i=hi+1;i<rows.length;i++) {
    const row=rows[i];if(row.every(v=>!v.trim()))continue;
    if(row.length!==header.length)throw new Error(`${i+1}行目の列数が一致しません。`);
    const numeric=[row[0],...indexes.map(j=>row[j])].map(v=>v.trim()===''?NaN:Number(v));
    if(!numeric.every(Number.isFinite))throw new Error(`${i+1}行目の時刻・加速度に欠損または非数値があります。補間せず解析を停止しました。`);
    time.push(numeric[0]);names.forEach((n,j)=>signals[n].push(numeric[j+1]));
  }
  if(time.length<3)throw new Error('時系列が3点未満です。');
  const steps=time.slice(1).map((t,i)=>t-time[i]);
  if(steps.some(d=>d<=0))throw new Error('時刻は重複のない昇順にしてください。');
  const sorted=[...steps].sort((a,b)=>a-b), mid=Math.floor(sorted.length/2);
  const dt=sorted.length%2?sorted[mid]:(sorted[mid-1]+sorted[mid])/2;
  if(steps.some(d=>Math.abs(d-dt)>1e-7+0.01*dt))throw new Error('時刻間隔が不均一です。サンプルの欠落・同期を確認してください。');
  const metadata=Object.fromEntries(rows.slice(0,heads[0]).filter(r=>r.length===2));
  if(metadata['記録周波数(Hz)']!==undefined && (!Number.isFinite(Number(metadata['記録周波数(Hz)']))||Math.abs(Number(metadata['記録周波数(Hz)'])-1/dt)>.01/dt))throw new Error('記録周波数と時刻列が一致しません。');
  const units=rows.slice(0,hi).filter(r=>r[0]==='単位'&&r.length===header.length).at(-1);
  if(units&&names.some(n=>!['mm/s/s','mm/s^2','mm/s²'].includes(units[header.indexOf(n)])))throw new Error('モード・頭部加速度の単位はmm/s²である必要があります。');
  return {name,time,signals,metadata,fs:1/dt,min:time[0],max:time.at(-1)+dt};
}

export function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
export function demean(a){const m=mean(a);return a.map(v=>v-m);}
export function rms(a){return Math.sqrt(mean(a.map(v=>v*v)));}
export function metrics(ankle,hip) {
  if(ankle.length!==hip.length||ankle.length<3||![...ankle,...hip].every(Number.isFinite))throw new Error('有効な同じ長さの波形が3点以上必要です。');
  const a=demean(ankle),h=demean(hip),ra=rms(a),rh=rms(h),total=ra+rh;
  return {ankle_rms:ra,hip_rms:rh,ankle_percent:total>1e-12?100*ra/total:null,
    hip_percent:total>1e-12?100*rh/total:null,rxi:ra>1e-12?rh/ra:null,
    cxi:Math.min(ra,rh)>1e-12?Math.max(-1,Math.min(1,mean(a.map((v,i)=>v*h[i]))/(ra*rh))):null};
}

export function analyze(recording,start=5,end=35) {
  if(!Number.isFinite(start)||!Number.isFinite(end)||start>=end)throw new Error('開始・終了に数値を入力し、開始 < 終了にしてください。');
  if(start<recording.min-1e-7||end>recording.max+1e-7)throw new Error(`解析区間を記録範囲（${recording.min.toFixed(2)}〜${recording.max.toFixed(2)}秒）内に設定してください。`);
  const ids=[];recording.time.forEach((t,i)=>{if(t>=start-1e-9&&t<end-1e-9)ids.push(i);});
  if(ids.length<3)throw new Error('選択区間に3点以上のデータが必要です。');
  const axes={},warnings=[];
  for(const [axis,[ac,hc,head]] of Object.entries(AXES)) {
    const a=ids.map(i=>recording.signals[ac][i]), h=ids.map(i=>recording.signals[hc][i]);
    const m=metrics(a,h);
    if(m.ankle_percent===null)warnings.push(`${axis}: 両モードの変動がゼロのため、割合は未定義です。`);
    if(m.rxi===null)warnings.push(`${axis}: 足関節RMSがゼロのため、RXIは未定義です。`);
    if(recording.signals[head]) {
      const hd=demean(ids.map(i=>recording.signals[head][i])), da=demean(a),dh=demean(h);
      m.head_rms=rms(hd);m.reconstruction_rmse=rms(hd.map((v,i)=>da[i]+dh[i]-v));
      m.reconstruction_relative_percent=m.head_rms>1e-12?100*m.reconstruction_rmse/m.head_rms:null;
      m.consistent=m.reconstruction_rmse<=Math.max(1e-5,.01*m.head_rms);
      if(!m.consistent)warnings.push(`${axis}: モード和と頭部加速度にずれがあります（平均除去後RMSE ${m.reconstruction_rmse.toFixed(3)} mm/s²）。割合はCSVのモード列から算出しています。頭部加速度への寄与を解釈する際は、機器側の計算定義を確認してください。`);
    } else {m.consistent=null;warnings.push(`${axis}: 頭部加速度列がないため、モード和の整合性は未確認です。`);}
    axes[axis]=m;
  }
  return {source:recording.name,reference:REFERENCE,start_s:start,end_s:end,samples:ids.length,
    sampling_hz:recording.fs,input:'exported_processed_modes',interval:'start inclusive, end exclusive',
    percentage_definition:'100 * RMS(mode) / (RMS(ankle) + RMS(hip)); demeaned',
    additional_filter:false,axes,warnings,ids};
}

export function resultCSV(result) {
  const q=s=>'"'+String(s??'').replaceAll('"','""')+'"';
  // Prefix untrusted filename cells to prevent spreadsheet formula execution.
  const filename=/^[=+\-@\t\r]/.test(result.source)?"'"+result.source:result.source;
  const rows=[['source','start_s','end_s','samples','axis','ankle_rms_mm_s2','hip_rms_mm_s2','ankle_rms_share_percent','hip_rms_share_percent','RXI','CXI','reconstruction_rmse_mm_s2','consistency','notes']];
  for(const [axis,m]of Object.entries(result.axes))rows.push([filename,result.start_s,result.end_s,result.samples,axis,...['ankle_rms','hip_rms','ankle_percent','hip_percent','rxi','cxi','reconstruction_rmse'].map(k=>m[k]),m.consistent===null?'unverified':m.consistent?'within_tolerance':'check_required',result.warnings.join(' / ')]);
  return '\uFEFF'+rows.map(row=>row.map(q).join(',')).join('\r\n');
}
