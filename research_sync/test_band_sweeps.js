const fs=require('fs'),path=require('path'),assert=require('assert');
const w=require('./waveform.js');
for(let band=1;band<=3;band++)for(const direction of ['up','down']){
 const id=`study_sweep100_band${band}_${direction}`;
 const c=JSON.parse(fs.readFileSync(path.join(__dirname,'study_profiles',id+'.json')));
 assert.deepStrictEqual(c.analysis_interval_s,[10,95]);
 assert.strictEqual(c.duration_s,100);
 assert.strictEqual(c.hold_final_position,true);
 assert((direction==='up')===(c.start_frequency_hz<c.end_frequency_hz));
 const k=Math.log(c.end_frequency_hz/c.start_frequency_hz)/85;
 let worstV=0,worstA=0;
 for(let t=.013;t<99.99;t+=.137){
  const h=.0001,l=w.evaluate(t-h,c),r=w.evaluate(t+h,c),q=w.evaluate(t,c);
  for(const value of Object.values(q))assert(Number.isFinite(value));
  assert(Math.abs(q.v)<=.1+1e-12);
  worstV=Math.max(worstV,Math.abs((r.x-l.x)/(2*h)-q.v));
  worstA=Math.max(worstA,Math.abs((r.v-l.v)/(2*h)-q.a));
  if(t>10&&t<95){
   const phase=2*Math.PI*c.start_frequency_hz*10+
    2*Math.PI*c.start_frequency_hz*Math.expm1(k*(t-10))/k;
   assert(Math.abs(q.v-.1*Math.sin(phase))<1e-12);
  }
 }
 assert(worstV<2e-7,`${id} velocity derivative ${worstV}`);
 assert(worstA<1e-7,`${id} acceleration derivative ${worstA}`);
 for(const t of [10,95]){
  const l=w.evaluate(t-1e-7,c),r=w.evaluate(t+1e-7,c);
  assert(Math.abs(l.v-r.v)<1e-7);assert(Math.abs(l.a-r.a)<1e-7);
 }
 assert.deepStrictEqual(w.evaluate(-1,c),{x:0,v:0,a:0});
 assert.deepStrictEqual(w.evaluate(0,c),{x:0,v:0,a:0});
 const end=w.evaluate(100,c);
 assert.strictEqual(end.v,0);assert.strictEqual(end.a,0);
 assert.deepStrictEqual(w.evaluate(101,c),end);
 const before=w.evaluate(100-1e-4,c);
 assert(Math.abs(before.v)<1e-10);assert(Math.abs(before.a)<1e-7);
 // Independent midpoint integral of the analytic velocity checks displacement.
 const dt=.002;let integral=0;
 for(let i=0;i<50000;i++)integral+=w.evaluate((i+.5)*dt,c).v*dt;
 assert(Math.abs(integral-end.x)<1e-7,`${id} integral mismatch`);
 const invalid={...c,onset_s:0};assert.throws(()=>w.evaluate(0,invalid),/Invalid/);
 console.log(`PASS ${id}: continuous x/v/a, derivatives, velocity envelope, integral, final hold`);
}
