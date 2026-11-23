"""Overlay four 20s cycles of stimulus and unconverted force/moment."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import detrend
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'cycle_repeatability'

def correlations(b):
    c=np.corrcoef(b)
    return c[np.triu_indices(4,1)]

def main():
    OUT.mkdir(exist_ok=True);rows=[];pairrows=[];saved=[]
    fig,axes=plt.subplots(3,4,figsize=(17,10),sharex=True)
    rng=np.random.default_rng(20261124)
    # Keep clear of boundary extrapolation; 50 Hz interpolation grid, no smoothing.
    phase=np.arange(0.,20.,.02)
    for col,folder in enumerate(sorted(ROOT.glob('20261121_*_90s_multicosine_*'))):
        stem=folder.name;q=pd.read_csv(folder/f'{stem}_quest.csv');d=pd.read_csv(folder/f'{stem}_combalance.csv')
        qt=q.pc_mapped_time_s.to_numpy();qe=q.stimulus_elapsed_s.to_numpy();ft=d.pc_estimated_sample_time_s.to_numpy()
        endpts=np.interp([10,90],qe,qt)
        assert ft.min()<=endpts[0] and ft.max()>=endpts[1]-.1
        starts=[];targets=[]
        for k in range(4):
            elapsed=10+20*k+phase;pc=np.interp(elapsed,qe,qt);targets.append(pc)
        signals={
            'VR_velocity':np.array([np.interp(pc,qt,q.ap_velocity_m_s) for pc in targets]),
            'sensor_Fy':np.array([np.interp(pc,ft,d.sensor_Fy_N) for pc in targets]),
            'sensor_Mx':np.array([np.interp(pc,ft,d.sensor_Mx_Nm) for pc in targets])}
        title=f'{q.subject.iloc[0]} trial {int(q.trial_number.iloc[0])}'
        local,la=plt.subplots(3,1,figsize=(11,8),sharex=True)
        for row,(name,raw) in enumerate(signals.items()):
            b=raw-raw.mean(axis=1,keepdims=True);avg=b.mean(axis=0)
            pair=correlations(b);dtpair=correlations(detrend(b,axis=1))
            freq=np.fft.rfftfreq(len(phase),.02)
            spec=np.fft.rfft(b,axis=1)
            band=(freq>=.04999)&(freq<=.70001)
            bb=np.fft.irfft(spec*band, n=len(phase),axis=1)
            bandpair=correlations(bb)
            loo=np.array([np.corrcoef(b[k],np.mean(np.delete(b,k,axis=0),axis=0))[0,1] for k in range(4)])
            obs=pair.mean()
            # Shift each whole cycle independently: retains waveform/autocorrelation.
            # Exploratory circular-shift reference, not a definitive hypothesis test.
            null=[]
            if name!='VR_velocity':
                for _ in range(2000):
                    shifted=np.array([np.roll(v,int(s)) for v,s in zip(b,rng.integers(0,len(phase),4))])
                    null.append(correlations(shifted).mean())
                p=(1+sum(v>=obs for v in null))/(len(null)+1)
            else:p=float('nan')
            rows.append(dict(trial=stem,channel=name,mean_pair_correlation=obs,min_pair_correlation=min(pair),max_pair_correlation=max(pair),mean_pair_correlation_linear_detrended=dtpair.mean(),mean_pair_correlation_0p05_0p70Hz=bandpair.mean(),mean_leave_one_cycle_out_correlation=loo.mean(),cycle_average_rms=float(np.sqrt(np.mean(avg**2))),single_cycle_rms=float(np.sqrt(np.mean(b**2))),exploratory_shift_p=p,cycles=4))
            n=0
            for i in range(4):
                for j in range(i+1,4):pairrows.append(dict(trial=stem,channel=name,cycle_a=i+1,cycle_b=j+1,correlation=pair[n]));n+=1
            for i,t in enumerate(phase):
                saved.append(dict(trial=stem,channel=name,cycle_time_s=t,cycle_average=avg[i],**{f'cycle{k+1}':b[k,i] for k in range(4)}))
            unit={'VR_velocity':'m/s','sensor_Fy':'N','sensor_Mx':'Nm'}[name]
            for ax in [axes[row,col],la[row]]:
                for k in range(4):ax.plot(phase,b[k],alpha=.5,lw=.7,label=f'Cycle {k+1}')
                ax.plot(phase,avg,color='black',lw=1.7,label='Cycle average')
                ax.grid(alpha=.25);ax.set_xlim(0,20)
            axes[row,col].text(.02,.96,f'Mean pair r = {obs:.2f}',transform=axes[row,col].transAxes,va='top',fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
            if col==0:axes[row,col].set_ylabel(f'{name} [{unit}]\ncycle mean removed')
            la[row].set_ylabel(f'{name} [{unit}]\ncycle mean removed')
        axes[0,col].set_title(title);axes[-1,col].set_xlabel('Time within 20s cycle [s]')
        la[0].legend(ncol=5,fontsize=8,loc='upper right');la[-1].set_xlabel('Time within 20s cycle [s]')
        local.suptitle(title+': raw stimulus/force/moment, four cycles (no COM conversion)')
        local.tight_layout();local.savefig(OUT/f'{stem}_overlay.png',dpi=170);plt.close(local)
    axes[0,0].legend(fontsize=7,ncol=3,loc='lower right')
    fig.suptitle('Four 20s cycles: stimulus and raw sensor channels; offsets removed separately per cycle')
    fig.tight_layout();fig.savefig(OUT/'all_trials_overlay.png',dpi=170);plt.close(fig)
    pd.DataFrame(rows).to_csv(OUT/'repeatability_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(pairrows).to_csv(OUT/'pair_correlations.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(saved).to_csv(OUT/'aligned_cycles.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
