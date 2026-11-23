// Execute scripts extracted from the delivered HTML, not separately injected modules.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync(__dirname+'/two-link-lab.html','utf8');
assert(!/<script\b[^>]*\bsrc\s*=/i.test(html),'all executable code is embedded');
const elements={};
for(const tag of html.matchAll(/<([a-z][\w-]*)\b([^>]*)>/gi)){
 const id=tag[2].match(/\bid="([^"]+)"/);if(!id)continue;
 assert(!elements[id[1]],'unique element IDs');
 const attrs=tag[2],value=attrs.match(/\bvalue="([^"]*)"/);
 elements[id[1]]={value:value?value[1]:'',checked:/\bchecked\b/.test(attrs),style:{},hidden:false,textContent:'',setAttribute(){},classList:{toggle(){}},setPointerCapture(){},getBoundingClientRect:()=>({width:800,height:470})};
 if(tag[1]==='select')elements[id[1]].value=html.slice(tag.index).match(/<option\s+value="([^"]+)"/)[1];
}
const ctx=new Proxy({}, {get:(_,name)=>name==='arc'? (...a)=>{assert(a.every(Number.isFinite));assert(a[2]>=0)}:()=>{},set:()=>true});
for(const el of Object.values(elements))el.getContext=()=>ctx;
function boot(withModule){
 let next,time=1;const listeners={};
 const env={document:{getElementById:id=>elements[id]||null,querySelectorAll:()=>[],addEventListener(){}},requestAnimationFrame:f=>next=f,addEventListener:(type,f)=>listeners[type]=f};
 if(withModule)env.module={exports:{}}; // Some embedded viewers define this global.
 vm.createContext(env);
 for(const script of html.matchAll(/<script>([\s\S]*?)<\/script>/gi))vm.runInContext(script[1],env);
 const frames=n=>{for(let i=0;i<n;i++){time+=1000/60;next(time)}};
 frames(3);assert.equal(elements['startup-status'].hidden,true,'first rendered frame clears loading message');
 assert.equal(elements.profile.value,withModule?'reference':'older');
 assert.equal(typeof elements['world-reference'].onpointerdown,'function','left model accepts shared input');
 elements.auto.onclick();frames(60);assert.notEqual(elements.floor.textContent,'0.0 cm');
 assert.equal(elements.floor.textContent,elements['floor-reference'].textContent,'both models receive identical floor motion');
 for(const profile of ['older','gainOnly','delayOnly','reference']){elements.profile.value=profile;elements.profile.onchange();frames(120);assert(elements.comparison.textContent.includes('RMS'));assert.equal(elements.floor.textContent,elements['floor-reference'].textContent);if(profile==='reference')assert.equal(elements.hip.textContent,elements['hip-reference'].textContent,'same profile produces matching angles')}
 listeners.error({error:new Error('test diagnostic')});assert(!elements['startup-status'].hidden);assert(elements['startup-status'].textContent.includes('test diagnostic'));
 elements.profile.value='reference';
}
boot(false);boot(true);
assert.equal(fs.readFileSync(__dirname+'/index.html','utf8'),html,'both entry points contain the same release');
console.log('PASS: delivered HTML has no external scripts; real markup and inline scripts initialize, animate, switch all profiles, and display errors, with/without a module global. Browser rendering not tested.');
