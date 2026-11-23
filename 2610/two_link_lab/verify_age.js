const assert=require('node:assert/strict'),model=require('./physics.js');
const base={floor:0,scene:0,control:true,weight:.5},dt=1/240;
// The original reference trajectory must remain bit-for-bit identical.
let plain=Array(8).fill(0),sim=model.createSimulation();
for(let i=0;i<2400;i++){const input={...base,floor:.01*Math.sin(i/120)};plain=model.step(plain,input);assert.deepEqual(sim.step(input),plain)}
// Scene perturbations cannot affect body angles before the additional delay elapses.
sim=model.createSimulation('delayOnly');let s;
for(let i=0;i<5;i++){s=sim.step({...base,scene:.02});assert.equal(s[0],0);assert.equal(s[1],0)}
s=sim.step({...base,scene:.02});assert(Math.hypot(s[2],s[3])>0,'feedback begins after 24 ms');
function run(id,dt,seconds=20,forcing=t=>({...base,scene:.02*Math.sin(2*Math.PI*1.2*t)})){
 const sim=model.createSimulation(id);let state;
 for(let t=0;t<seconds-dt/2;t+=dt){state=sim.step(forcing(t),dt);assert(state.every(Number.isFinite));assert(Math.max(...state.slice(0,2).map(Math.abs))<20*Math.PI/180)}
 return state;
}
for(const id of Object.keys(model.profiles)){
 const profile=model.profiles[id];
 const held=run(id,dt,30,()=>({...base,scene:.02}));
 for(let j=0;j<2;j++)assert(Math.abs(held[j]-model.params.V[j]*.5*.02/(profile.kScale*model.params.K[j]))<1e-5,'correct static gain with delayed feedback');
 const a=run(id,dt),b=run(id,dt/2),c=run(id,dt/4);
 const err=(u,v)=>Math.max(...u.map((x,i)=>Math.abs(x-v[i])));
 assert(err(a,b)<.001,'time-step convergence');assert(err(b,c)<err(a,b),'refinement reduces numerical error');
 const blind=run(id,dt,3,t=>({...base,scene:.02*Math.sin(t*5),weight:0}));
 assert(blind.slice(0,6).every(v=>v===0),'visual weight zero works for every profile');
 console.log(id,'PASS, step refinement error:',err(b,c).toExponential(3));
}
// With active control off, the age-profile controller changes must have no effect.
const free=id=>run(id,dt,.4,t=>({...base,floor:.002*Math.sin(t*3),control:false}));
assert.deepEqual(free('reference'),free('older'));
assert.throws(()=>model.createSimulation('unknown'));
assert.throws(()=>model.createSimulation('older').step(base,.1));
console.log('PASS: reference regression, delay causality, static gain, numerical refinement, finite responses, sensory decoupling, control-off equivalence.');
