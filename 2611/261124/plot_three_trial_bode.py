"""Measured harmonic FRFs for the three Nov 24 equal-velocity trials."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SOURCE = Path(__file__).resolve().parents[2]/'research_sync'
PRIOR = SOURCE/'results/equal_velocity_20261124'
OUT = Path(__file__).resolve().parent/'three_trial_bode_sign_checked'
FREQ = np.array([.05,.10,.20,.30,.45,.60,.70,.90,1.20,1.50])
EXCITED=FREQ.copy()
SUBJECTS = ['seguchi','otsubo','takato']
NAMES = ['瀬口','大坪','たかと']
ORDERS = ['変位','速度','加速度']
START,END=10.,90.

def amplitude(t,y):
    mask = np.isfinite(t)&np.isfinite(y)&(t>=START)&(t<END)
    t,y = np.asarray(t)[mask],np.asarray(y)[mask]
    assert len(t)>1000 and t.min()<START+.03 and t.max()>END-.05
    # Fit all excited frequencies and the nuisance offset/drift jointly.
    design=np.column_stack([np.ones(len(t)),t-t.mean()]+
        [v for f in EXCITED for v in (np.cos(2*np.pi*f*t),np.sin(2*np.pi*f*t))])
    beta=np.linalg.lstsq(design,y,rcond=None)[0]
    return (beta[2::2]-1j*beta[3::2])[np.isin(EXCITED,FREQ)]

def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'Yu Gothic','font.size':9,'axes.unicode_minus':False})
    results={}; rows=[]; s=2j*np.pi*FREQ
    for subject in SUBJECTS:
        folder=SOURCE/'recordings/saved'/f'20261124_{subject}_100s_multicosine10_equal_velocity_001'
        q=pd.read_csv(folder/f'{folder.name}_quest.csv')
        head=pd.read_csv(PRIOR/'gravity_correction'/f'{subject}_gravity_corrected.csv')
        b=pd.read_csv(folder/f'{folder.name}_combalance.csv')
        meta=json.loads((folder/f'{folder.name}_metadata.json').read_text(encoding='utf-8-sig'))
        mass=meta['weight_estimate']['body_mass_kg'];height=meta['participant']['height_cm']/100
        mb=.978*mass;lb=.513*height;jb=.0473*mass*height**2
        coefficient=jb/(mb*lb)+.039*height+lb
        origin_below=-float(meta['device_parameters']['FpAz0'])
        # User confirms same facing as the right-heel-drop reference.
        # OEM quiet/AP swaps X/Y force and moment channels without negation.
        # Match standing_analysis.py: Ny_contact=Mx-origin_below*Fy.
        com_position=(b.sensor_Mx_Nm-(origin_below+coefficient)*b.sensor_Fy_N)/(mb*9.80665)
        elapsed=np.interp(b.pc_estimated_sample_time_s,q.pc_mapped_time_s,q.stimulus_elapsed_s,left=np.nan,right=np.nan)
        inputs=[amplitude(q.stimulus_elapsed_s,q[c]) for c in
                ['ap_displacement_m','ap_velocity_m_s','ap_acceleration_m_s2']]
        ha=amplitude(head.stimulus_elapsed_s,head.corrected_accel_tau30_m_s2)
        cx=amplitude(elapsed,com_position)
        for kind,outputs in [('head',[ha/s**2,ha/s,ha]),('com',[cx,s*cx,s**2*cx])]:
            for i,u in enumerate(inputs):
                assert np.all(abs(u)>1e-9)
                for j,y in enumerate(outputs):
                    h=y/u
                    assert np.isfinite(h).all()
                    results[kind,subject,i,j]=h
                    for f,z in zip(FREQ,h):
                        rows.append(dict(output_source=kind,subject=subject,input_order=i,
                            output_order=j,frequency_hz=f,real=z.real,imag=z.imag,
                            gain=abs(z),gain_db=20*np.log10(abs(z)),phase_deg=np.angle(z,deg=True)))
    pd.DataFrame(rows).to_csv(OUT/'frequency_responses.csv',index=False)
    for subject,name in zip(SUBJECTS,NAMES):
        fig,axes=plt.subplots(6,3,figsize=(15,17),sharex=True,layout='constrained')
        for j in range(3):
            for i in range(3):
                ag,ap=axes[2*j,i],axes[2*j+1,i]
                for kind,label,color in [('head','頭部IMU','#0072B2'),('com','COM推定','#D55E00')]:
                    h=results[kind,subject,i,j]
                    ag.semilogx(FREQ,20*np.log10(abs(h)),'o-',ms=3,label=label,color=color)
                    ap.semilogx(FREQ,np.angle(h,deg=True),'o-',ms=3,color=color)
                ag.set_title(f'入力：刺激{ORDERS[i]} → 出力：{ORDERS[j]}')
                unit=['m','m/s','m/s²']
                ag.set_ylabel(f'ゲイン [dB re 1 ({unit[j]})/({unit[i]})]')
                ap.set_ylabel('位相 [deg]');ap.set_ylim(-190,190);ap.set_yticks([-180,-90,0,90,180])
                for ax in (ag,ap):ax.grid(True,which='both',alpha=.25)
        axes[0,0].legend(ncol=2)
        for ax in axes[-1]:ax.set_xlabel('周波数 [Hz]')
        fig.suptitle(f'2026/11/24 {name} | 頭部IMUとCOM推定の比較\n{START:g}–{END:g} s | 頭部：暫定重力補正 τ=30 s・COM：既存男性身体係数 | 位相：機器間遅延未校正',fontsize=13)
        fig.savefig(OUT/f'{subject}_head_com_bode_3x3.png',dpi=150)
        plt.close(fig)
    for kind,title in [('head','頭部IMU（暫定重力補正 τ=30 s）'),('com','COM推定位置（既存の男性身体係数）')]:
        fig,axes=plt.subplots(6,3,figsize=(15,17),sharex=True,layout='constrained')
        for j in range(3):
            for i in range(3):
                ag,ap=axes[2*j,i],axes[2*j+1,i]
                for subject,label in zip(SUBJECTS,NAMES):
                    h=results[kind,subject,i,j]
                    ag.semilogx(FREQ,20*np.log10(abs(h)),'o-',ms=3,label=label)
                    ap.semilogx(FREQ,np.angle(h,deg=True),'o-',ms=3)
                ag.set_title(f'入力：刺激{ORDERS[i]} → 出力：{ORDERS[j]}')
                unit=['m','m/s','m/s²']
                ag.set_ylabel(f'ゲイン [dB re 1 ({unit[j]})/({unit[i]})]')
                ap.set_ylabel('位相 [deg]');ap.set_ylim(-190,190);ap.set_yticks([-180,-90,0,90,180])
                for ax in (ag,ap):ax.grid(True,which='both',alpha=.25)
        axes[0,0].legend(ncol=3)
        for ax in axes[-1]:ax.set_xlabel('周波数 [Hz]')
        fig.suptitle(f'2026/11/24 等速度マルチサイン 3試行 | {title}\n解析区間 {START:g}–{END:g} s・励起10周波数の同時最小二乗 | 線は測定点の接続',fontsize=14)
        fig.savefig(OUT/f'{kind}_bode_3x3.png',dpi=150)
        plt.close(fig)
    (OUT/'method.json').write_text(json.dumps(dict(
        subjects=SUBJECTS,window_s=[START,END],frequencies_hz=FREQ.tolist(),
        input='Three logged stimulus channels; joint least squares includes all 10 excited frequencies as nuisance terms, report only selected frequencies. Intercept and trend included.',
        head='Existing calibrated, provisional gravity-corrected acceleration tau=30s; velocity=A/(iw), displacement=A/(iw)^2 at each excited frequency. No 1Hz position filter applied.',
        com='Recomputed unfiltered position, male coefficients: Rx=Fy, Ny_contact=Mx-0.0123*Fy; X=(Ny_contact-c*Fy)/(mb*g). Same facing as heel-drop reference confirmed by user. Velocity=iwX, acceleration=(iw)^2X.',
        phase='Wrapped [-180,180]. PC clock mapping, device delay uncalibrated.',
        limitations='Head gravity compensation provisional. OEM AP mapping and origin correction follow standing_analysis.py; same physical facing confirmed by user, horizontal-force known-load calibration still absent. Mass/anthropometry provisional. One trial per participant; no confidence intervals or coherence inferred.'
    ),ensure_ascii=False,indent=2),encoding='utf-8')
    print(OUT)
    print(f'{len(rows)} frequency response points; all finite')

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--start',type=float,default=10);parser.add_argument('--end',type=float,default=90);parser.add_argument('--output',type=Path);parser.add_argument('--max-frequency',type=float,default=1.5)
    args=parser.parse_args();START,END=args.start,args.end
    assert END>START and (END-START)>=20 and abs((END-START)/20-round((END-START)/20))<1e-8
    if args.output:OUT=args.output.resolve()
    FREQ=EXCITED[EXCITED<=args.max_frequency+1e-10]
    main()
