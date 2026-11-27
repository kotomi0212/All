const fs=require('fs'),path=require('path');
const root=__dirname,source=path.resolve(root,'../../research_sync');
const w=require(path.join(source,'waveform.js'));
const out=path.join(root,'sweep100');fs.mkdirSync(out,{recursive:true});
const summary=[];
for(let band=1;band<=3;band++)for(const direction of ['up','down']){
 const id=`study_sweep100_band${band}_${direction}`;
 const c=JSON.parse(fs.readFileSync(path.join(source,'study_profiles',id+'.json')));
 const rows=['time_s,displacement_m,velocity_m_s,acceleration_m_s2,frequency_hz,segment'];
 const peak={x:0,v:0,a:0};
 for(let i=0;i<=10000;i++){
  const t=i/100,q=w.evaluate(t,c),u=Math.max(0,Math.min(85,t-10));
  const f=c.start_frequency_hz*Math.exp(Math.log(c.end_frequency_hz/c.start_frequency_hz)*u/85);
  const segment=t<10?'onset':t<95?'sweep':t<100?'offset':'hold';
  rows.push([t,q.x,q.v,q.a,f,segment].join(','));
 }
 // Peak checks use the actual 1 ms integration grid, not the coarser plot CSV.
 for(let i=0;i<=100000;i++){
  const q=w.evaluate(i/1000,c);
  for(const key of ['x','v','a'])peak[key]=Math.max(peak[key],Math.abs(q[key]));
 }
 fs.writeFileSync(path.join(out,id+'.csv'),rows.join('\n'));
 fs.writeFileSync(path.join(out,id+'.json'),JSON.stringify(c,null,2));
 summary.push({condition:id,start_frequency_hz:c.start_frequency_hz,end_frequency_hz:c.end_frequency_hz,
  peak_abs_displacement_m:peak.x,peak_abs_velocity_m_s:peak.v,peak_abs_acceleration_m_s2:peak.a,
  final_position_m:w.evaluate(100,c).x});
}
fs.writeFileSync(path.join(out,'summary.json'),JSON.stringify(summary,null,2));
console.log(JSON.stringify(summary,null,2));
