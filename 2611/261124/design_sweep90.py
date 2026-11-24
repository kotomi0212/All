"""Full 90s constant-acceleration downsweep, with no extra trial phases."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from design_visual_stimuli import A,metrics

OUT=Path(__file__).resolve().parent/'stimulus_design'/'sweep90'
T,F0,F1=90.,.70,.05

def build(dt):
    t=np.arange(round(T/dt)+1)*dt;k=np.log(F1/F0)/T
    f=F0*np.exp(k*t);theta=2*np.pi*F0*np.expm1(k*t)/k
    c=np.trapezoid(np.cos(theta),t);s=np.trapezoid(np.sin(theta),t)
    phase=np.arctan2(c,s)
    # Phase is selected before integration; no subsequent drift correction.
    a=A*np.cos(theta+phase);v=cumulative_trapezoid(a,t,initial=0);x=cumulative_trapezoid(v,t,initial=0)
    return t,x,v,a,f,float(phase)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    t,x,v,a,f,p=build(.001);fine=build(.0005)
    np.testing.assert_allclose(x,fine[1][::2],atol=2e-6,rtol=0)
    np.testing.assert_allclose(v,fine[2][::2],atol=2e-7,rtol=0)
    assert abs(v[-1])<1e-10 and np.isclose(f[0],F0) and np.isclose(f[-1],F1)
    k=np.log(F1/F0)/T;theta=2*np.pi*F0*np.expm1(k*t)/k+p
    np.testing.assert_allclose(a,A*np.cos(theta),atol=1e-12)
    cfg=dict(status='design_only_not_installed_in_VR',duration_s=T,analysis_interval_s=[0,T],start_frequency_hz=F0,end_frequency_hz=F1,acceleration_amplitude_m_s2=A,phase_rad=p,frequency_formula='f(t)=0.70*exp(log(0.05/0.70)*t/90)',acceleration_formula='a(t)=A*cos(2*pi*0.70*expm1(k*t)/k+phase), k=log(0.05/0.70)/90',onset_s=0,offset_s=0,return_s=0,phase_selection='zero net acceleration integral, v(0)=v(90)=0; no drift subtraction',after_trial='Hold the final room position; do not reset position to zero at 90s',maxima=metrics(x,v,a),initial_state=dict(x=float(x[0]),v=float(v[0]),a=float(a[0])),final_state=dict(x=float(x[-1]),v=float(v[-1]),a=float(a[-1])),total_oscillations=float((F1-F0)/k),integration_max_position_difference_m=float(max(abs(x-fine[1][::2]))),note='Constant acceleration cosine envelope over the entire 90s, no ramps. Position/velocity join rest continuously; acceleration generally steps at start/end. Low-frequency response is transient, not stationary FRF.')
    pd.DataFrame(dict(time_s=t,ap_displacement_m=x,ap_velocity_m_s=v,ap_acceleration_m_s2=a,instantaneous_frequency_hz=f)).to_csv(OUT/'downsweep90.csv',index=False)
    (OUT/'downsweep90.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(4,1,figsize=(11,9),sharex=True)
    for ax,y,label in zip(axes,[x,v,a,f],['Displacement [m]','Velocity [m/s]','Acceleration [m/s2]','Frequency [Hz]']):
        ax.plot(t[::10],y[::10],lw=.9);ax.set(ylabel=label,xlim=(0,90));ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time [s]');fig.suptitle('Full 90s downsweep: 0.70 -> 0.05 Hz; constant acceleration amplitude; no extra phases')
    fig.tight_layout();fig.savefig(OUT/'downsweep90.png',dpi=170);plt.close(fig)
    (OUT/'notes.txt').write_text('90秒全体を0.70→0.05Hzの対数下降スイープとする。立ち上げ・減衰停止・原点復帰の別区間はない。加速度のコサイン振幅は全区間で約0.0167m/s^2。初期位相を選び、初速・終速ともゼロ。終了後は最終位置で保持し、原点へ瞬間リセットしない。位置・速度は停止状態につながるが、加速度は開始・終了で段差を持つ。開始過渡と低周波側の滞在時間不足を含むため、全90秒を定常FRFとして扱わない。高周波を戻して当初指定範囲を網羅する一方、各周波数で定常応答を測れることは保証しない。振幅正規化・後付けのドリフト除去なし。既存VRアプリは未変更。\n',encoding='utf-8')
    print(json.dumps(cfg,indent=2))

if __name__=='__main__':main()
