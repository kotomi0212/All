const M=require('./model.js'),C=M.complex;
function diff(f){const r=M.response(f);let d=(C.phase(r.theta[0])-C.phase(r.theta[1]))*180/Math.PI;while(d>180)d-=360;while(d<=-180)d+=360;return d}
let ok=true;for(const f of [.024,.1,.2,.8,1.2,1.5,2.936]){const r=M.response(f),d=diff(f);console.log(`${f.toFixed(3)} Hz  |theta1|=${C.abs(r.theta[0]).toExponential(3)} rad/u  |theta2|=${C.abs(r.theta[1]).toExponential(3)} rad/u  phase(theta1-theta2)=${d.toFixed(1)} deg`);if(!Number.isFinite(d))ok=false}
const low=Math.abs(diff(.024)),high=Math.abs(diff(2.936));if(low>30||high<100){console.error(`Qualitative validation failed: low=${low.toFixed(1)}°, high=${high.toFixed(1)}°`);ok=false}else console.log(`PASS: low frequency nearly in-phase (${low.toFixed(1)}°); high-frequency phase separation (${high.toFixed(1)}°).`);
process.exitCode=ok?0:1;
