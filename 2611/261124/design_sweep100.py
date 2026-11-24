"""100s visual acceleration downsweep draft; no VR application changes."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from design_visual_stimuli import A,smooth,metrics

OUT=Path(__file__).resolve().parent/'stimulus_design'/'sweep100'
START_F,END_F=.30,.05
ONSET,SWEEP,OFFSET,RETURN=5.,70.,15.,10.

def build(dt):
    T=ONSET+SWEEP+OFFSET;t=np.arange(round(T/dt)+1)*dt
    u=t-ONSET;k=np.log(END_F/START_F)/SWEEP;core=np.clip(u,0,SWEEP)
    phi=2*np.pi*START_F*np.expm1(k*core)/k
    phi=np.where(u<0,2*np.pi*START_F*u,phi+2*np.pi*END_F*np.maximum(u-SWEEP,0))
    f=np.where(u<0,START_F,np.where(u<=SWEEP,START_F*np.exp(k*core),END_F))
    e=smooth(t/ONSET)*smooth((T-t)/OFFSET)
    c=np.trapezoid(e*np.cos(phi),t);s=np.trapezoid(e*np.sin(phi),t);phase=np.arctan2(c,s)
    a=A*e*np.cos(phi+phase);v=cumulative_trapezoid(a,t,initial=0);x=cumulative_trapezoid(v,t,initial=0)
    rt=np.arange(1,round(RETURN/dt)+1)*dt;z=rt/RETURN;end=x[-1]
    rx=end*(1-smooth(z));rv=-end*(30*z*z-60*z**3+30*z**4)/RETURN
    ra=-end*(60*z-180*z*z+120*z**3)/RETURN**2
    return np.r_[t,T+rt],np.r_[x,rx],np.r_[v,rv],np.r_[a,ra],np.r_[f,np.full(len(rt),np.nan)],phase

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    t,x,v,a,f,phase=build(.001);fine=build(.0005)
    np.testing.assert_allclose(x,fine[1][::2],atol=2e-6,rtol=0)
    np.testing.assert_allclose(v,fine[2][::2],atol=2e-7,rtol=0)
    assert abs(v[90000])<1e-10
    np.testing.assert_allclose([x[0],v[0],a[0],x[-1],v[-1],a[-1]],0,atol=1e-10)
    core=(t>=5)&(t<=75);k=np.log(END_F/START_F)/SWEEP
    theta=2*np.pi*START_F*np.expm1(k*(t[core]-ONSET))/k+phase
    np.testing.assert_allclose(a[core],A*np.cos(theta),atol=1e-12)
    config=dict(status='design_only_not_installed_in_VR',duration_s=100,start_frequency_hz=START_F,end_frequency_hz=END_F,sweep_duration_s=SWEEP,analysis_interval_s=[5,75],onset_interval_s=[0,5],offset_interval_s=[75,90],return_interval_s=[90,100],acceleration_amplitude_m_s2=A,phase_rad=float(phase),frequency_formula='f(t)=0.30*exp(log(0.05/0.30)*(t-5)/70), 5<=t<=75',acceleration_formula='a(t)=A*cos(2*pi*0.30*expm1(k*(t-5))/k+phase), k=log(0.05/0.30)/70',integration='1ms cumulative trapezoid; 0.5ms convergence checked',maxima=metrics(x,v,a),max_abs_acceleration_analysis_m_s2=float(max(abs(a[core]))),return_start_position_m=float(x[90000]),cycles_in_sweep=float((END_F-START_F)/k),note='Amplitude is constant only in 5..75s. Onset/offset and return are excluded. Short sweep is not quasi-static at all frequencies. No post-hoc drift subtraction or amplitude renormalization.')
    pd.DataFrame(dict(time_s=t,ap_displacement_m=x,ap_velocity_m_s=v,ap_acceleration_m_s2=a,instantaneous_frequency_hz=f,analysis=core)).to_csv(OUT/'downsweep100.csv',index=False)
    (OUT/'downsweep100.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(4,1,figsize=(11,9),sharex=True)
    for ax,y,label in zip(axes,[x,v,a,f],['Displacement [m]','Velocity [m/s]','Acceleration [m/s2]','Frequency [Hz]']):
        ax.plot(t[::10],y[::10],lw=.9);ax.axvspan(5,75,color='green',alpha=.10);ax.axvspan(90,100,color='gray',alpha=.12)
        ax.set(ylabel=label,xlim=(0,100));ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time [s]');fig.suptitle(f'100s downsweep: 0.30 -> 0.05 Hz; A={A:.7f} m/s2; green = analysis')
    fig.tight_layout();fig.savefig(OUT/'downsweep100.png',dpi=170);plt.close(fig)
    print(json.dumps(config,indent=2))

if __name__=='__main__':main()
