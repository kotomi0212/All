const assert=require('node:assert/strict');
const model=require('./physics.js'),p=model.params,zero=()=>Array(8).fill(0);
const base={floor:0,scene:0,control:true,weight:.5};
function run(seconds,signal,state=zero(),dt=1/240){for(let t=0;t<seconds-dt/2;t+=dt)state=model.step(state,signal(t),dt);return state}
assert(p.M[0][0]>0&&p.M[0][0]*p.M[1][1]-p.M[0][1]**2>0,'positive kinetic energy');
assert.deepEqual(run(5,()=>base),zero(),'upright equilibrium');
let s=zero();s[0]=.02;s[1]=-.01;
const stable=run(20,()=>base,s);
assert(Math.max(...stable.slice(0,4).map(Math.abs))<1e-5,'recover from initial perturbation');
const falling=run(.5,()=>({...base,control:false}),s);
assert(Math.hypot(...falling.slice(0,2))>Math.hypot(...s.slice(0,2)),'uncontrolled inverted body diverges');
const first=model.derivative(zero(),{...base,floor:.02});
const generalized=p.M.map(r=>r[0]*first[2]+r[1]*first[3]);
generalized.forEach((v,i)=>assert(Math.abs(v+p.H[i]*first[5])<1e-10,'base acceleration inertial balance'));
assert(first[2]<0,'positive floor acceleration initially tilts leg backward');
const invisible=run(10,t=>({...base,weight:0,scene:.03*Math.sin(8*t)}));
assert(invisible.slice(0,6).every(v=>v===0),'zero visual weight decouples scene from body');
const target=run(20,()=>({...base,floor:.04}));
assert(Math.abs(target[4]-.04)<1e-6,'floor smoothly reaches held cursor position');
assert(Math.max(...target.slice(0,4).map(Math.abs))<.001,'body recovers after floor stops');
const signal=t=>({...base,floor:.02*Math.sin(2*Math.PI*1.2*t)});
const fine=run(4,signal,zero(),1/480),coarse=run(4,signal);
assert(Math.max(...fine.map((v,i)=>Math.abs(v-coarse[i])))<.002,'time-step convergence with sampled command');
for(const f of [.2,1.2]){
 let a=0,b=0,ab=0,state=zero();
 for(let i=0;i<14400;i++){
  state=model.step(state,{...base,scene:2*Math.PI/180*Math.sin(2*Math.PI*f*i/240)});
  assert(state.every(Number.isFinite));
  if(i>7200){a+=state[0]**2;b+=state[1]**2;ab+=state[0]*state[1]}
 }
 const correlation=ab/Math.sqrt(a*b);
 console.log(`Illustrative visual response ${f} Hz: leg/HAT correlation ${correlation.toFixed(3)}`);
 assert(f===.2?correlation>.9:correlation<-.5,'documented qualitative example');
}
console.log('PASS: equilibrium, recovery, instability, floor inertia, visual decoupling, held floor, convergence, frequency example.');
