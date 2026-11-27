(function(root){
  const apDownsweepCache=new WeakMap();
  const bandSweepCache=new WeakMap();

  function evaluateBandSweep(t,c){
    let data=bandSweepCache.get(c);
    if(!data){
      if(c.duration_s!==100||c.onset_s!==10||c.sweep_s!==85||c.offset_s!==5||
         c.sample_step_s!==.001||c.velocity_amplitude_m_s!==.1||
         c.initial_position_m!==0||c.initial_phase_rad!==0||c.hold_final_position!==true||
         !Number.isFinite(c.start_frequency_hz)||!Number.isFinite(c.end_frequency_hz)||
         Math.min(c.start_frequency_hz,c.end_frequency_hz)<.05-1e-12||
         Math.max(c.start_frequency_hz,c.end_frequency_hz)>.5+1e-12||
         c.start_frequency_hz===c.end_frequency_hz)throw new Error('Invalid band-sweep configuration');
      const start=c.onset_s,end=start+c.sweep_s,T=c.duration_s;
      const f0=c.start_frequency_hz,f1=c.end_frequency_hz;
      const k=Math.log(f1/f0)/c.sweep_s;
      const phase0=2*Math.PI*f0*start;
      const phaseEnd=phase0+2*Math.PI*f0*Math.expm1(k*c.sweep_s)/k;
      function state(time){
        let phase,f;
        if(time<start){phase=2*Math.PI*f0*time;f=f0;}
        else if(time<=end){const u=time-start;phase=phase0+2*Math.PI*f0*Math.expm1(k*u)/k;f=f0*Math.exp(k*u);}
        else{phase=phaseEnd+2*Math.PI*f1*(time-end);f=f1;}
        let e=1,ed=0;
        if(time<=0||time>=T){e=0;}
        else if(time<start||time>end){
          const ramp=time<start?start:c.offset_s;
          const u=time<start?time/ramp:(T-time)/ramp;
          e=u*u*u*(10+u*(-15+6*u));
          ed=30*u*u*(1-u)*(1-u)/ramp*(time<start?1:-1);
        }
        return {v:c.velocity_amplitude_m_s*e*Math.sin(phase),
          a:c.velocity_amplitude_m_s*(ed*Math.sin(phase)+e*2*Math.PI*f*Math.cos(phase))};
      }
      const dt=c.sample_step_s,n=Math.round(T/dt),positions=new Float64Array(n+1);
      let prev=state(0).v;
      for(let i=1;i<=n;i++){
        const next=state(i*dt).v;
        positions[i]=positions[i-1]+(prev+next)*dt/2;prev=next;
      }
      data={positions,state};bandSweepCache.set(c,data);
    }
    const {positions,state}=data,T=c.duration_s,dt=c.sample_step_s,last=positions.length-1;
    if(t<=0)return{x:positions[0],v:0,a:0};
    if(t>=T)return{x:positions[last],v:0,a:0};
    const q=state(t),i=Math.min(last-1,Math.floor(t/dt));
    const x=positions[i]+(state(i*dt).v+q.v)*(t-i*dt)/2;
    return{x,v:q.v,a:q.a};
  }

  function evaluateApDownsweep(t,c){
    let data=apDownsweepCache.get(c);
    const dt=c.sample_step_s;
    if(!data){
      const k=Math.log(c.end_frequency_hz/c.start_frequency_hz)/c.duration_s;
      const omega0=2*Math.PI*c.start_frequency_hz;
      const phaseAt=time=>c.initial_phase_rad+omega0*Math.expm1(k*time)/k;
      const intervals=Math.round(c.duration_s/dt);
      const positions=new Float64Array(intervals+1);
      positions[0]=c.initial_position_m;
      let previousVelocity=c.velocity_amplitude_m_s*Math.sin(c.initial_phase_rad);
      for(let i=1;i<=intervals;i++){
        const velocity=c.velocity_amplitude_m_s*Math.sin(phaseAt(i*dt));
        positions[i]=positions[i-1]+(previousVelocity+velocity)*dt/2;
        previousVelocity=velocity;
      }
      data={positions,phaseAt,k,omega0};
      apDownsweepCache.set(c,data);
    }
    const {positions,phaseAt,k,omega0}=data;
    const last=positions.length-1;
    if(t<0)return{x:positions[0],v:0,a:0};
    if(t>c.duration_s)return{x:positions[last],v:0,a:0};

    const phase=phaseAt(t);
    const v=c.velocity_amplitude_m_s*Math.sin(phase);
    const a=c.velocity_amplitude_m_s*omega0*Math.exp(k*t)*Math.cos(phase);
    if(t===c.duration_s)return{x:positions[last],v,a};
    const i=Math.min(last-1,Math.floor(t/dt));
    const sampleTime=i*dt;
    const sampleVelocity=c.velocity_amplitude_m_s*Math.sin(phaseAt(sampleTime));
    // Integrate the fractional last step when a display frame falls between samples.
    const x=positions[i]+(sampleVelocity+v)*(t-sampleTime)/2;
    return{x,v,a};
  }

  function evaluate(t,c){
    if(c.waveform==='ap_velocity_log_band_sweep')return evaluateBandSweep(t,c);
    if(c.waveform==='ap_velocity_log_downsweep')return evaluateApDownsweep(t,c);
    if(c.waveform==='sampled_centered_sweep'){
      const T=c.duration_s,h=c.sample_step_s,s=c.samples_xva;
      if(t<0)return{x:s[0][0],v:0,a:0};
      if(t>T)return{x:s[s.length-1][0],v:0,a:0};
      const i=Math.min(s.length-2,Math.floor(t/h)),u=(t-i*h)/h;
      const [x0,v0,a0]=s[i],[x1,v1,a1]=s[i+1];
      // Quintic Hermite interpolation: x, v, a are derivatives of the same curve.
      const c0=x0,c1=h*v0,c2=h*h*a0/2;
      const dx=x1-c0-c1-c2,dv=h*v1-c1-2*c2,da=h*h*a1-2*c2;
      const c3=10*dx-4*dv+da/2,c4=-15*dx+7*dv-da,c5=6*dx-3*dv+da/2;
      return{x:c0+u*(c1+u*(c2+u*(c3+u*(c4+u*c5)))),
        v:(c1+u*(2*c2+u*(3*c3+u*(4*c4+u*5*c5))))/h,
        a:(2*c2+u*(6*c3+u*(12*c4+u*20*c5)))/(h*h)};
    }
    if(c.waveform==='pilot_acceleration_multicosine'){
      const T=c.duration_s,R=c.ramp_s;if(t<=0||t>=T)return{x:0,v:0,a:0};
      let x=0,v=0,a=0;
      c.frequencies_hz.forEach((f,k)=>{const w=2*Math.PI*f,z=w*t+c.phases_rad[k],amp=c.acceleration_components_m_s2?.[k]??c.acceleration_component_m_s2;x-=amp*Math.cos(z)/w**2;v+=amp*Math.sin(z)/w;a+=amp*Math.cos(z);});
      const u=Math.max(0,Math.min(1,Math.min(t,T-t)/R));
      const e=10*u**3-15*u**4+6*u**5,d=u>0&&u<1?(30*u**2-60*u**3+30*u**4)/R*(t<T/2?1:-1):0;
      const dd=u>0&&u<1?(60*u-180*u**2+120*u**3)/R**2:0;
      return{x:e*x,v:e*v+d*x,a:e*a+2*d*v+dd*x};
    }
    if(c.waveform==='fixation_single_sine')return (root.FixationStimulus||(typeof require!=='undefined'?require('./fixation.js'):null)).state(t);
    if(c.schema===2&&c.waveform==='endpoint_stationary_multicosine'){
      const T=c.duration_s,R=c.stationary_start_s;
      if(t<0||t>T)return{x:0,v:0,a:0};
      let x=0,v=0,a=0;
      c.frequencies_hz.forEach((f,k)=>{
        const w=2*Math.PI*f,theta=c.endpoint_phases_rad[k],z=w*(t-T)+theta;
        x+=(Math.cos(theta)-Math.cos(z))/w**2;v+=Math.sin(z)/w;a+=Math.cos(z);
      });
      const u=Math.max(0,Math.min(1,t/R)),e=10*u**3-15*u**4+6*u**5;
      const d=t>0&&t<R?(30*u**2-60*u**3+30*u**4)/R:0;
      const dd=t>0&&t<R?(60*u-180*u**2+120*u**3)/R**2:0;
      const A=c.acceleration_component_m_s2;
      return{x:A*e*x,v:A*(e*v+d*x),a:A*(e*a+2*d*v+dd*x)};
    }
    let x=0,v=0,a=0;
    c.frequencies_hz.forEach((f,k)=>{const w=2*Math.PI*f,z=w*t+c.phases_rad[k];x-=Math.cos(z)/w**2;v+=Math.sin(z)/w;a+=Math.cos(z);});
    const A=c.acceleration_component_m_s2,u=Math.max(0,Math.min(1,Math.min(t,90-t)/5));
    const e=10*u**3-15*u**4+6*u**5;
    const d=(u>0&&u<1)?(30*u**2-60*u**3+30*u**4)*(t<45?1:-1)/5:0;
    const dd=(u>0&&u<1)?(60*u-180*u**2+120*u**3)/25:0;
    return {x:A*e*x,v:A*(e*v+d*x),a:A*(e*a+2*d*v+dd*x)};
  }
  root.ResearchWaveform={evaluate};
  if(typeof module!=='undefined')module.exports={evaluate};
})(typeof window!=='undefined'?window:globalThis);
