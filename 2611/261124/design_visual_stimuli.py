"""Draft visual-motion stimuli; exports only, does not change the VR application."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

OUT=Path(__file__).resolve().parent/'stimulus_design'
A=.01666464946520325
F=np.array([.05,.10,.20,.30,.45,.60,.70,.90,1.20,1.50])

def smooth(u):
    u=np.clip(u,0,1)
    return 10*u**3-15*u**4+6*u**5

def envelope(t,T,R):
    u=np.clip(t/R,0,1);v=np.clip((T-t)/R,0,1)
    def parts(z):return smooth(z),(30*z*z-60*z**3+30*z**4)/R,(60*z-180*z*z+120*z**3)/R**2
    e,d,dd=parts(u);h,dh,ddh=parts(v)
    return e*h,d*h-e*dh,dd*h-2*d*dh+e*ddh

def multi(t,phases):
    w=2*np.pi*F;theta=t[:,None]*w+phases
    x=-A*np.sum(np.cos(theta)/w**2,axis=1)
    v=A*np.sum(np.sin(theta)/w,axis=1);a=A*np.cos(theta).sum(axis=1)
    e,d,dd=envelope(t,100,10)
    return e*x,e*v+d*x,e*a+2*d*v+dd*x

def sweep(dt):
    # 10s onset at 0.70Hz; 240s logarithmic downsweep; 20s offset at 0.05Hz.
    R,D,E=10.,240.,20.;T=R+D+E
    t=np.arange(round(T/dt)+1)*dt;u=t-R;k=np.log(.05/.70)/D
    core=np.clip(u,0,D)
    phase=2*np.pi*.70*np.expm1(k*core)/k
    phase=np.where(u<0,2*np.pi*.70*u,phase+2*np.pi*.05*np.maximum(u-D,0))
    freq=np.where(u<0,.70,np.where(u<=D,.70*np.exp(k*core),.05))
    e=smooth(t/R)*smooth((T-t)/E)
    C=np.trapezoid(e*np.cos(phase),t);S=np.trapezoid(e*np.sin(phase),t)
    phi=np.arctan2(C,S) # exact discrete integral of acceleration is zero
    acc=A*e*np.cos(phase+phi)
    vel=cumulative_trapezoid(acc,t,initial=0);pos=cumulative_trapezoid(vel,t,initial=0)
    # Return to neutral after the analysis, using a 20s quintic move.
    rt=np.arange(1,round(20/dt)+1)*dt;z=rt/20;end=pos[-1]
    rx=end*(1-smooth(z));rv=-end*(30*z*z-60*z**3+30*z**4)/20
    ra=-end*(60*z-180*z*z+120*z**3)/20**2
    return np.r_[t,T+rt],np.r_[pos,rx],np.r_[vel,rv],np.r_[acc,ra],np.r_[freq,np.full(len(rt),np.nan)],phi

def metrics(x,v,a):
    return dict(max_abs_displacement_m=float(max(abs(x))),max_abs_velocity_m_s=float(max(abs(v))),max_abs_acceleration_m_s2=float(max(abs(a))))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261124)
    # Reproducible phase selection; common component amplitude is never rescaled.
    tg=np.arange(0,100.0001,.025);candidates=rng.uniform(-np.pi,np.pi,(512,len(F)))
    def score(p):
        ys=multi(tg,p)
        return max(max(abs(y))/limit for y,limit in zip(ys,[.30,.10,.20]))
    phases=min(candidates,key=score)
    t=np.arange(100001)/1000;x,v,a=multi(t,phases)
    np.testing.assert_allclose([x[0],v[0],a[0],x[-1],v[-1],a[-1]],0,atol=1e-12)
    core=(t>=10)&(t<90);fc=np.fft.rfftfreq(core.sum(),.001);ac=2*abs(np.fft.rfft(a[core]))/core.sum()
    lines=ac[np.rint(F*80).astype(int)]
    np.testing.assert_allclose(lines,A,rtol=1e-10,atol=1e-12)
    # 20s periodicity in the pure interval, with exact requested component amplitudes.
    for yy in multi(np.arange(10,30,.01),phases):assert np.isfinite(yy).all()
    for first,second in zip(multi(np.arange(10,30,.01),phases),multi(np.arange(30,50,.01),phases)):
        np.testing.assert_allclose(first,second,atol=1e-11)
    mconfig=dict(status='design_only_not_installed_in_VR',waveform='ten_component_acceleration_multicosine',duration_s=100,analysis_interval_s=[10,90],ramp_s=10,period_s=20,analysis_cycles=4,frequencies_hz=F.tolist(),phases_rad=phases.tolist(),acceleration_component_m_s2=A,normalization='none',maxima=metrics(x,v,a),formula='x0=-sum(A*cos(2*pi*f*t+phi)/(2*pi*f)^2); x=e*x0; v=e*v0+e_dot*x0; a=e*a0+2*e_dot*v0+e_ddot*x0',note='Equal component acceleration only in the 10..90s core; startup/shutdown add sidebands. No claim of stronger SNR than the old stimulus.')
    pd.DataFrame(dict(time_s=t,ap_displacement_m=x,ap_velocity_m_s=v,ap_acceleration_m_s2=a,analysis=core)).to_csv(OUT/'multicosine10.csv',index=False)
    pd.DataFrame(dict(frequency_hz=F,acceleration_amplitude_m_s2=lines,phase_rad=phases)).to_csv(OUT/'multicosine10_components.csv',index=False)
    (OUT/'multicosine10.json').write_text(json.dumps(mconfig,indent=2),encoding='utf-8')
    st,sx,sv,sa,sf,phi=sweep(.001)
    fine=sweep(.0005)
    np.testing.assert_allclose(sx,fine[1][::2],atol=2e-6,rtol=0)
    np.testing.assert_allclose(sv,fine[2][::2],atol=2e-7,rtol=0)
    assert abs(sv[270000])<1e-10
    np.testing.assert_allclose([sx[0],sv[0],sa[0],sx[-1],sv[-1],sa[-1]],0,atol=1e-10)
    sc=(st>=10)&(st<=250)
    # Check exact acceleration generation rather than treating x=-A/w(t)^2*cos(phi) as exact.
    k=np.log(.05/.70)/240;theta=2*np.pi*.70*np.expm1(k*(st[sc]-10))/k+phi
    np.testing.assert_allclose(sa[sc],A*np.cos(theta),atol=1e-12)
    sconfig=dict(status='design_only_not_installed_in_VR',waveform='constant_acceleration_log_downsweep',duration_s=290,analysis_interval_s=[10,250],sweep_duration_s=240,start_frequency_hz=.70,end_frequency_hz=.05,onset_s=10,offset_s=20,return_to_neutral_s=20,acceleration_amplitude_m_s2=A,phase_rad=float(phi),frequency_formula='f=0.70*exp(log(0.05/0.70)*(t-10)/240) for 10<=t<=250',acceleration_formula='a=A*cos(phase), phase=2*pi*0.70*expm1(k*(t-10))/k+phi in core; acceleration envelope outside core',integration='cumulative trapezoid, 1ms; validated against 0.5ms',phase_selection='zero net acceleration impulse over 0..270s; no linear drift subtraction',maxima=metrics(sx,sv,sa),position_before_return_m=float(sx[270000]),sweep_cycles=float((.05-.70)/k),note='Not periodic. Reaching zero speed at 270s leaves a position offset; a separate quintic return ends at 290s. Analyze only 10..250s. No automatic amplitude normalization.')
    pd.DataFrame(dict(time_s=st,ap_displacement_m=sx,ap_velocity_m_s=sv,ap_acceleration_m_s2=sa,instantaneous_frequency_hz=sf,analysis=sc)).to_csv(OUT/'downsweep.csv',index=False)
    (OUT/'downsweep.json').write_text(json.dumps(sconfig,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(4,2,figsize=(14,11))
    for col,(tt,xx,vv,aa,title,interval) in enumerate([(st,sx,sv,sa,'Log downsweep 0.70 -> 0.05 Hz',[10,250]),(t,x,v,a,'10-component multicosine',[10,90])]):
        for row,(yy,label) in enumerate(zip([xx,vv,aa],['Displacement [m]','Velocity [m/s]','Acceleration [m/s2]'])):
            axes[row,col].plot(tt[::20],yy[::20],lw=.7);axes[row,col].axvspan(*interval,color='green',alpha=.08)
            axes[row,col].set(xlabel='Time [s]',ylabel=label);axes[row,col].grid(alpha=.2)
        axes[0,col].set_title(title)
    axes[3,0].plot(st,sf);axes[3,0].set(xlabel='Time [s]',ylabel='Frequency [Hz]');axes[3,0].axvspan(10,250,color='green',alpha=.08)
    axes[3,1].stem(F,lines);axes[3,1].set(xlabel='Frequency [Hz]',ylabel='Component acceleration [m/s2]',ylim=(0,A*1.2))
    fig.suptitle(f'Draft visual stimulus designs; A = {A:.7f} m/s2; green = analysis interval')
    fig.tight_layout();fig.savefig(OUT/'design_overview.png',dpi=170);plt.close(fig)
    (OUT/'verification.json').write_text(json.dumps(dict(multicosine=mconfig,sweep=sconfig,max_relative_multicosine_amplitude_error=float(max(abs(lines/A-1))),sweep_displacement_max_integration_difference_m=float(max(abs(sx-fine[1][::2])))),indent=2),encoding='utf-8')
    print(json.dumps(dict(multicosine=mconfig,sweep=sconfig),indent=2))

if __name__=='__main__':main()
