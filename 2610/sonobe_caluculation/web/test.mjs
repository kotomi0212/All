import assert from 'node:assert/strict';
import fs from 'node:fs';
import {parseCSV,readRecording,analyze,metrics,resultCSV} from './dist/analysis.mjs';
assert.deepEqual(parseCSV('a,"b,b","c""d"\r\n1,2,3\r\n'),[['a','b,b','c"d'],['1','2','3']]);
assert.throws(()=>parseCSV('a,"b'),/引用符/);
const a=[-1,0,1],m=metrics(a.map(v=>v+100),a.map(v=>-2*v+200));
assert.ok(Math.abs(m.ankle_percent-100/3)<1e-10);assert.equal(m.rxi,2);assert.ok(Math.abs(m.cxi+1)<1e-12);
assert.equal(metrics([0,0,0],[0,0,0]).ankle_percent,null);
assert.equal(metrics([0,0,0],a).rxi,null);
const csv=['記録周波数(Hz),100',',xix1,xix2,xiy1,xiy2,axhd,ayhd','単位,mm/s/s,mm/s/s,mm/s/s,mm/s/s,mm/s/s,mm/s/s','MIN,0,0,0,0,0,0',',xix1,xix2,xiy1,xiy2,axhd,ayhd'];
for(let i=0;i<4000;i++){let t=i/100,a=Math.sin(t*2*Math.PI),h=2*Math.cos(t*2*Math.PI);csv.push([t,a,h,a,h,a+h,a+h].join(','));}
const r=readRecording(csv.join('\n'),'=test.csv'),res=analyze(r);
assert.equal(res.samples,3000);assert.ok(Math.abs(res.axes.AP.rxi-2)<1e-10);assert.equal(res.warnings.length,0);
assert.throws(()=>analyze(r,-1,35));assert.throws(()=>analyze(r,5,41));assert.throws(()=>analyze(r,NaN,35));
assert.throws(()=>readRecording('time,Fx,AccX1\n0,1,1'),/処理済み/);
assert.throws(()=>readRecording(csv.join('\n').replace('0,0,2,0,2,2,2','0,nan,2,0,2,2,2')),/欠損/);
assert.ok(resultCSV(res).includes("'=test.csv"));
r.signals.ayhd=r.signals.ayhd.map(v=>v*2);assert.ok(analyze(r).warnings.some(w=>w.startsWith('ML:')));
if(process.argv[2]){
 const raw=fs.readFileSync(process.argv[2]),rec=readRecording(new TextDecoder('shift_jis',{fatal:true}).decode(raw),'Test_00010.csv'),out=analyze(rec);
 const expected={AP:[40.8651074000588,47.586840784219,1.164485885681302,-.46512935450040394],ML:[56.7968305934706,73.9066267884559,1.3012456155775751,-.5264012370483062]};
 for(const [axis,values]of Object.entries(expected))['ankle_rms','hip_rms','rxi','cxi'].forEach((k,i)=>assert.ok(Math.abs(out.axes[axis][k]-values[i])<1e-10));
 assert.equal(out.axes.AP.consistent,true);assert.equal(out.axes.ML.consistent,false);
 console.log('Sample: all 8 exported indices match; ML consistency warning verified.');
}
for(const name of ['index.html','app.mjs','analysis.mjs','style.css'])assert.ok(fs.existsSync(new URL('./dist/'+name,import.meta.url)));
const client=fs.readFileSync(new URL('./dist/app.mjs',import.meta.url),'utf8');
assert.ok(!/\b(fetch|XMLHttpRequest|WebSocket|sendBeacon|localStorage|indexedDB)\b/.test(client));
console.log('Web analysis tests passed: CSV, offsets, ratios, zero, interval, malformed data, warnings, exports.');
