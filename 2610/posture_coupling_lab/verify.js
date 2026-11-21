global.window=global;require('./model.js');const M=global.PostureModel,C=M.complex;let ok=true;
for(const inertia of [true,false])for(const control of [false,true])for(const crossScale of [0,1,2])for(let i=0;i<80;i++){const f=.1+1.9*i/79,r=M.response(f,{inertia,control,crossScale});for(const z of [...r.q,r.hip])if(!Number.isFinite(C.abs(z)))ok=false}
const coupled=M.response(.8,{inertia:true,control:false}),decoupled=M.response(.8,{inertia:false,control:false});
const A=M.matrices({inertia:true,control:false}),B=M.matrices({inertia:true,control:true});
if(A.Kjoint[0][1]!==0||A.Kjoint[1][0]!==0||B.Kjoint[0][1]===0||B.Kjoint[1][0]===0)ok=false;
const slow=M.response(.2,{control:true}),mechanicsOff=M.response(.8,{control:false,inertia:false}),controlStrong=M.response(.8,{control:true,crossScale:1.5}),hipPreset=M.response(1.4,{control:true,crossScale:1.3,passiveScale:.9});
let slowPhase=Math.abs((C.phase(slow.q[0])-C.phase(slow.q[1]))*180/Math.PI);while(slowPhase>180)slowPhase=Math.abs(slowPhase-360);
const mechanicsChange=Math.abs(C.abs(coupled.hip)-C.abs(mechanicsOff.hip))/C.abs(coupled.hip),hipRatio=C.abs(hipPreset.hip)/C.abs(hipPreset.q[0]);
if(slowPhase>10||mechanicsChange<.15||hipRatio<1.2||Math.abs(C.abs(controlStrong.q[0])-C.abs(coupled.q[0]))<1e-4)ok=false;
console.log('0.8 Hz Model A hip gain, inertia ON:',C.abs(coupled.hip).toFixed(6));console.log('0.8 Hz Model A hip gain, inertia OFF:',C.abs(decoupled.hip).toFixed(6));console.log('KJ,A off-diagonal:',A.Kjoint[0][1],A.Kjoint[1][0]);console.log('KJ,B off-diagonal:',B.Kjoint[0][1],B.Kjoint[1][0]);console.log('preset checks: slow phase diff',slowPhase.toFixed(1),'deg; inertia effect',(mechanicsChange*100).toFixed(1),'%; hip/ankle gain',hipRatio.toFixed(2));console.log(ok?'PASS: responses, decomposition, and recommended presets satisfy their stated behavior.':'FAIL');process.exitCode=ok?0:1;
