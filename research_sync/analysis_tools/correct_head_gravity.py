"""Provisional 6-axis gravity compensation; no true-position claim."""
import csv,json,argparse
from pathlib import Path
import numpy as np
from scipy.signal import butter,sosfiltfilt,detrend
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

G=9.80665
def read(p,names):
    with p.open(encoding='utf-8-sig',newline='') as f:r=list(csv.DictReader(f))
    return {k:np.array([float(x[k]) for x in r]) for k in names}
def compensate(acc,gyro,dt,pre,tau):
    bias=np.median(gyro[pre],axis=0)
    up=acc[pre].mean(axis=0);up=up/np.linalg.norm(up)
    ap=np.array([0.,0.,1.]);ap-=up*np.dot(ap,up);ap/=np.linalg.norm(ap)
    R=np.eye(3); output=[];tilt=[]
    for a,w in zip(acc,gyro-bias):
        R=R@Rotation.from_rotvec(w*dt).as_matrix()
        m=R@a
        if np.isfinite(tau):
            error=np.cross(m/np.linalg.norm(m),up)
            R=Rotation.from_rotvec(error*dt/tau).as_matrix()@R
        output.append(np.dot(R@a-G*up,ap))
        tilt.append(np.degrees(np.arcsin(np.clip(np.dot(R[:,2],up),-1,1))))
    return np.array(output),np.array(tilt),bias
def displacement(a,fs):
    centered=detrend(a)
    p=np.pad(centered,(len(a),len(a)),mode='reflect')
    f=np.fft.rfftfreq(len(p),1/fs);z=np.fft.rfft(p);band=(f>=.05)&(f<=4)
    x=np.zeros_like(z);x[band]=-z[band]/(2*np.pi*f[band])**2
    x=np.fft.irfft(x,n=len(p))[len(a):2*len(a)]
    return sosfiltfilt(butter(4,1,fs=fs,output='sos'),x)*1000

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trials',nargs='+',type=Path,help='Named recording directories; order defines anonymous A/B/C labels')
    parser.add_argument('--output',required=True,type=Path,help='New local output directory (must be git-ignored)')
    args=parser.parse_args()
    OUT=args.output.resolve()
    if OUT.exists():raise FileExistsError('Choose a new output directory')
    OUT.mkdir(parents=True)
    fig,axes=plt.subplots(len(args.trials),1,figsize=(15,3.3*len(args.trials)),sharex=True,layout='constrained',squeeze=False)
    audit=[]
    for index,(ax,folder) in enumerate(zip(axes.flat,args.trials)):
        s='Participant '+chr(65+index)
        stem=folder.name
        meta=json.loads((folder/f'{stem}_metadata.json').read_text(encoding='utf-8-sig'));params=meta['device_parameters']
        q=read(folder/f'{stem}_quest.csv',['pc_mapped_time_s','stimulus_elapsed_s','ap_displacement_m'])
        names=['pc_estimated_sample_time_s','sensor_Fy_N','sensor_Mx_Nm']+[f'external_imu_count_{i}' for i in range(6)]
        b=read(folder/f'{stem}_combalance.csv',names);ft=b[names[0]]
        fs=1/np.median(np.diff(ft));t=ft[0]+np.arange(int((ft[-1]-ft[0])*fs)+1)/fs
        if np.any(np.diff(ft)<=0) or np.max(np.diff(ft))>.1:raise ValueError('timestamps')
        elapsed=np.interp(t,q['pc_mapped_time_s'],q['stimulus_elapsed_s'],left=np.nan,right=np.nan)
        start=np.interp(0,q['stimulus_elapsed_s'],q['pc_mapped_time_s'])
        pre=(t>=start-2)&(t<start)
        if pre.sum()<100:raise ValueError('Insufficient initial baseline')
        raw=np.column_stack([np.interp(t,ft,b[f'external_imu_count_{i}']) for i in range(6)])
        ar=['G2','G4','G8','G16'].index(params['AccelRange']);gr=['D250','D500','D1000','D2000'].index(params['GyroRange'])
        ap=np.array(params['AccelPrf']['values']).reshape(2,3,4)[1,:,ar]
        az=np.array(params['AccelZero']['values']).reshape(2,3,4)[1,:,ar]
        gp=np.array(params['GyroPrf']['values']).reshape(2,3,4)[1,:,gr]
        gz=np.array(params['GyroZero']['values']).reshape(2,3,4)[1,:,gr]
        acc=(raw[:,:3]-az)*ap*G;gyro=np.deg2rad((raw[:,3:]-gz)*gp)
        ai=np.array(params['AccelItf']['values']).reshape(2,4,3,3)[1,ar]
        gi=np.array(params['GyroItf']['values']).reshape(2,4,3,3)[1,gr]
        if not np.allclose(ai,np.eye(3)) or not np.allclose(gi,np.eye(3)):raise ValueError('Nonidentity transform needs explicit OEM ordering audit')
        curves=[];accelerations=[];taus=[10.,30.,60.]
        valid=(elapsed>=10)&(elapsed<90)
        for tau in taus:
            a,tilt,bias=compensate(acc,gyro,1/fs,pre,tau)
            x=displacement(a,fs);x-=x[valid].mean();curves.append(x);accelerations.append(a)
        left,=ax.plot(q['stimulus_elapsed_s'],q['ap_displacement_m']*1000,':',color='gray',label='Visual stimulus (mm)')
        r=ax.twinx();r.fill_between(elapsed,np.min(curves,axis=0),np.max(curves,axis=0),color='#0072B2',alpha=.15,label='Head sensitivity: tau=10/30/60 s')
        line,=r.plot(elapsed,curves[1],color='#0072B2',lw=1,label='Gravity-corrected head estimate (tau=30 s)')
        mass=meta['weight_estimate']['body_mass_kg'];height=meta['participant']['height_cm']/100
        mb=.978*mass;lb=.513*height;jb=.0473*mass*height**2
        coefficient=jb/(mb*lb)+.039*height+lb
        cx=(-np.interp(t,ft,b['sensor_Mx_Nm'])-coefficient*np.interp(t,ft,b['sensor_Fy_N']))/(mb*G)*1000
        cx=sosfiltfilt(butter(4,1,fs=fs,output='sos'),cx)
        cx=cx[valid]-cx[valid].mean()
        com={'stimulus_elapsed_s':elapsed[valid],'sonobe_com_male_centered_mm':cx}
        cr=ax.twinx();cr.spines['right'].set_position(('axes',1.10))
        cl,=cr.plot(com['stimulus_elapsed_s'],com['sonobe_com_male_centered_mm'],color='#d62728',lw=1,label='Sonobe COM estimate')
        anonymous_label=s
        ax.set(title=anonymous_label,xlim=(0,100),ylabel='Visual (mm)');r.set_ylabel('Provisional head displacement (mm)',color='#0072B2');cr.set_ylabel('COM estimate (mm)',color='#d62728');ax.grid(alpha=.2)
        if index==0:ax.legend([left,line,cl],[a.get_label() for a in [left,line,cl]],fontsize=7,loc='upper right')
        ranges=[float(np.ptp(x[valid])) for x in curves]
        audit.append(dict(subject=s,taus_s=taus,head_peak_to_peak_mm_10_90=ranges,gyro_bias_deg_s=np.rad2deg(bias).tolist(),initial_accel_norm_g=float(np.linalg.norm(acc[pre].mean(axis=0))/G),head_sensor_Z_elevation_deg_range=[float(tilt[valid].min()),float(tilt[valid].max())]))
        with (OUT/f'{s}_gravity_corrected.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.writer(f);w.writerow(['stimulus_elapsed_s','head_mm_tau10','head_mm_tau30','head_mm_tau60','corrected_accel_tau30_m_s2']);w.writerows(zip(elapsed,*curves,accelerations[1]))
    axes[-1,0].set_xlabel('Stimulus elapsed time (s)')
    fig.suptitle('Provisional gyro + acceleration gravity correction | 1 Hz LP, below 0.05 Hz removed | blue shade: parameter sensitivity',fontsize=11)
    fig.savefig(OUT/'visual_head_com_gravity_corrected.png',dpi=150)
    report=dict(status='provisional compensation; not validated translation',method='6-axis complementary orientation, gyro propagation with accelerometer gravity alignment; tau=30 s, sensitivity 10/60 s. Initial AP = initial sensor Z projected onto horizontal. Gyro bias from 2 s before stimulus, assumed steady. OEM parameter scales and identity axis matrices. No fitted amplitude rescaling.',limitations=['Accelerometer tilt correction cannot distinguish sustained translation from gravity','Initial baseline not independently verified stationary','Channel axes inferred from installed OEM AdToEng; scale units supported by gravity test and gyro/gravity rotational agreement, no independent gyro fixture test','No magnetometer or independent heading; yaw and gyro bias drift remain','Hard lower frequency cutoff and reflected extension affect position','Smaller amplitude is not proof of accuracy'],results=audit)
    (OUT/'audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__ == "__main__":
    main()
