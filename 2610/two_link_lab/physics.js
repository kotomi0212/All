/* Small-angle sagittal two-link model. SI units; absolute angles from vertical.
   Anthropometry: Kiemel et al. 2008, Table 1. Controller is illustrative. */
(function(root){
 'use strict';
 const m=70.6,ma=m*.093,mb=m*.200,m2=m*.678,la=.420,lb=.378,l1=la+lb,l2=.511;
 const ha=la*.567,hb=lb*.567,h2=l2*.626,m1=ma+mb,h1=(ma*ha+mb*(la+hb))/m1;
 const J1=ma*(la*.643)**2+mb*(lb*.653)**2+mb*((la+hb)**2-hb**2),J2=m2*(l2*.798)**2;
 const M=[[J1+m2*l1*l1,m2*l1*h2],[m2*l1*h2,J2]],H=[m1*h1+m2*l1,m2*h2],G=H.map(v=>9.81*v);
 const K=[220,90],D=[70,12],V=[220,45],det=M[0][0]*M[1][1]-M[0][1]**2;
 // Wiesmeier et al. (2015), Results / Model parameters, Fig. 8:
 // P/mgh young=1.43, elderly=1.29; Td young=.153 s, elderly=.177 s.
 // Applying their whole-body gain ratio to both links and ONLY the delay difference
 // is a sensitivity experiment, not a fitted age-specific two-link controller.
 const ageEvidence=Object.freeze({url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC4477145/',youngP:1.43,olderP:1.29,youngDelay:.153,olderDelay:.177});
 const profiles=Object.freeze({
   reference:Object.freeze({label:'基準（従来の教育モデル）',kScale:1,extraDelay:0}),
   older:Object.freeze({label:'高齢者参考：P比＋遅延差',kScale:ageEvidence.olderP/ageEvidence.youngP,extraDelay:ageEvidence.olderDelay-ageEvidence.youngDelay}),
   gainOnly:Object.freeze({label:'要因分離：P比のみ',kScale:ageEvidence.olderP/ageEvidence.youngP,extraDelay:0}),
   delayOnly:Object.freeze({label:'要因分離：遅延差のみ',kScale:1,extraDelay:ageEvidence.olderDelay-ageEvidence.youngDelay})
 });
 function solve(v){return [(M[1][1]*v[0]-M[0][1]*v[1])/det,(-M[1][0]*v[0]+M[0][0]*v[1])/det]}
 // State = [leg angle, HAT angle, angular velocities, floor x, floor v, scene angle, scene velocity].
 function generalizedTorque(s,input,offset=0){
   const feedback=input.feedback?input.feedback(offset):s,ratio=input.kScale??1;
   // Artificial, instantaneous gravity compensation is retained in every profile.
   // Only the residual PD + visual command receives the extra feedback delay.
   return [0,1].map(i=>input.control ? -G[i]*s[i]-ratio*K[i]*feedback[i]+V[i]*input.weight*feedback[6]-D[i]*feedback[i+2] : 0);
 }
 function derivative(s,input,offset=0){
   const a=Math.max(-3,Math.min(3,35**2*(input.floor-s[4])-2*35*s[5]));
   const va=25**2*(input.scene-s[6])-2*25*s[7];
   const q=generalizedTorque(s,input,offset);
   const acc=solve(q.map((v,i)=>v+G[i]*s[i]-H[i]*a));
   return [s[2],s[3],acc[0],acc[1],s[5],a,s[7],va];
 }
 function step(s,input,dt=1/240){
   const k1=derivative(s,input),k2=derivative(s.map((v,i)=>v+dt*k1[i]/2),input,dt/2),k3=derivative(s.map((v,i)=>v+dt*k2[i]/2),input,dt/2),k4=derivative(s.map((v,i)=>v+dt*k3[i]),input,dt);
   return s.map((v,i)=>v+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6);
 }
 function torques(s,input){const Q=generalizedTorque(s,input);return [Q[0]+Q[1],Q[1]]}
 function createSimulation(profileId='reference',initial=Array(8).fill(0)){
   const profile=profiles[profileId];if(!profile)throw new Error('Unknown profile: '+profileId);
   let state=initial.slice(),time=0;const samples=[{t:0,s:state.slice()}];
   function at(t){
     if(t<=0)return initial;
     for(let i=1;i<samples.length;i++)if(samples[i].t>=t){const a=samples[i-1],b=samples[i],f=(t-a.t)/(b.t-a.t);return a.s.map((v,j)=>v+f*(b.s[j]-v))}
     throw new Error('Delay history does not cover requested time');
   }
   return {profile,step(input,dt=1/240){
     if(profile.extraDelay>0&&dt>profile.extraDelay)throw new Error('Step must not exceed feedback delay');
     const start=time,forcing={...input,kScale:profile.kScale,feedback:profile.extraDelay>0?offset=>at(start+offset-profile.extraDelay):null};
     state=step(state,forcing,dt);time+=dt;samples.push({t:time,s:state.slice()});
     while(samples.length>2&&samples[1].t<time-profile.extraDelay-dt)samples.shift();
     return state.slice();
   }};
 }
 const api={params:{m,m1,m2,l1,l2,h1,h2,J1,J2,M,H,G,K,D,V},profiles,ageEvidence,createSimulation,step,derivative,torques};
 root.TwoLink=api;
 if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(globalThis);
