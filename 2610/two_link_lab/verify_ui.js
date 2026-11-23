// Test application event wiring without a browser or rendering engine.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const ctx=new Proxy({}, {get:(_,name)=>name==='measureText'?()=>({width:10}):()=>{},set:()=>true});
const ids=['world-reference','selected-title','leg-reference','hat-reference','hip-reference','floor-reference','profile','profile-note','overlay','comparison','mode','position','frequency','amplitude','weight','control','pause','reset','manual','auto','world','trace','parameters','hint','status','leg','hat','hip','floor','position-value','frequency-value','amplitude-value','weight-value'];
const elements=Object.fromEntries(ids.map(id=>[id,{value:'0',checked:true,textContent:'',classList:{toggle(){}},getBoundingClientRect:()=>({width:800,height:id==='trace'?180:470}),getContext:()=>ctx,setPointerCapture(){}}]));
Object.assign(elements.mode,{value:'floor'});elements.profile.value='reference';elements.frequency.value='.2';elements.amplitude.value='2';elements.weight.value='.5';
let next;
const env={TwoLink:require('./physics.js'),devicePixelRatio:1,requestAnimationFrame:f=>next=f,document:{getElementById:id=>elements[id],querySelectorAll:()=>[],addEventListener(){}}};
vm.runInNewContext(fs.readFileSync(__dirname+'/app.js','utf8'),env);
let time=1;function frames(n){for(let i=0;i<n;i++){time+=1000/60;next(time)}}
frames(10);assert.equal(elements.leg.textContent,'0.00°');
elements.world.onpointerdown({clientX:400,pointerId:1});elements.world.onpointermove({clientX:500});elements.world.onpointerup();frames(120);
assert(+elements.floor.textContent.split(' ')[0]>0,'drag moves floor');
elements.pause.onclick();const held=elements.floor.textContent;frames(30);assert.equal(elements.floor.textContent,held,'pause freezes physics');
elements.reset.onclick();elements.mode.value='visual';elements.mode.onchange();elements.weight.value='0';elements.auto.onclick();frames(120);
assert.equal(elements.floor.textContent,'0.0 cm');assert.equal(elements.leg.textContent,'0.00°');assert.equal(elements.hat.textContent,'0.00°');
elements.weight.value='.5';frames(240);assert.notEqual(elements.hat.textContent,'0.00°');
elements.reset.onclick();frames(1);assert.equal(elements.floor.textContent,'0.0 cm');
for(const profile of ['older','gainOnly','delayOnly','reference']){
 elements.profile.value=profile;elements.profile.onchange();frames(1);
 assert.equal(elements.floor.textContent,'0.0 cm');
 frames(120);assert(!elements.status.textContent.includes('停止'));
 assert(elements.comparison.textContent.includes('RMS'));assert(elements.comparison.textContent.includes(env.TwoLink.profiles[profile].label));
}
console.log('PASS: application executes, drag, pause, reset, visual mode, visual weight, all profile switches, comparison readout. Rendering not tested.');

