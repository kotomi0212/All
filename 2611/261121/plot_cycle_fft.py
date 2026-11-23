"""One-sided FFT amplitudes of each 20s cycle, stimulus frequencies marked."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'cycle_repeatability'

def main():
    d=pd.read_csv(OUT/'aligned_cycles.csv'); rows=[]
    fig,axes=plt.subplots(3,4,figsize=(18,10),sharex=True)
    channels=['VR_velocity','sensor_Fy','sensor_Mx']
    labels=['VR velocity amplitude [m/s]','Fy amplitude [N]','Mx amplitude [Nm]']
    for col,trial in enumerate(d.trial.unique()):
        title=trial.replace('20261121_','').replace('_90s_multicosine_',' trial ')
        local,aa=plt.subplots(3,1,figsize=(12,9),sharex=True)
        for row,(channel,label) in enumerate(zip(channels,labels)):
            s=d[(d.trial==trial)&(d.channel==channel)].sort_values('cycle_time_s')
            dt=float(np.median(np.diff(s.cycle_time_s)));n=len(s)
            assert n==1000 and np.isclose(n*dt,20)
            f=np.fft.rfftfreq(n,dt)
            for k in range(1,5):
                y=s[f'cycle{k}'].to_numpy();y=y-y.mean()
                a=2*abs(np.fft.rfft(y))/n;a[0]/=2
                if n%2==0:a[-1]/=2
                keep=(f>0)&(f<=3.00001)
                for ax in [axes[row,col],aa[row]]:
                    ax.plot(f[keep],a[keep],'.-',ms=2.5,lw=.8,alpha=.85,label=f'Cycle {k}',zorder=3)
                for ff,v in zip(f[keep],a[keep]):rows.append(dict(trial=trial,channel=channel,cycle=k,frequency_hz=ff,one_sided_amplitude=v,stimulated=bool(ff<=1.500001)))
            for ax in [axes[row,col],aa[row]]:
                for ff in np.arange(.05,1.500001,.05):
                    ax.axvline(ff,color='#c43c39',alpha=.24,lw=.65,zorder=0)
                ax.set(yscale='log',xlim=(0,3));ax.grid(axis='y',alpha=.2)
            aa[row].set_ylabel(label)
            if col==0:axes[row,col].set_ylabel(label)
        axes[0,col].set_title(title);axes[-1,col].set_xlabel('Frequency [Hz]')
        aa[0].legend(ncol=4,fontsize=9);aa[-1].set_xlabel('Frequency [Hz]')
        local.suptitle(title+' - FFT per 20s cycle; red lines: stimulated frequencies')
        local.tight_layout();local.savefig(OUT/f'{trial}_fft.png',dpi=170);plt.close(local)
    axes[0,0].legend(ncol=2,fontsize=8)
    fig.suptitle('20s-cycle FFT amplitude (no averaging); red lines: 0.05 to 1.50 Hz in 0.05 Hz steps')
    fig.tight_layout();fig.savefig(OUT/'all_trials_fft.png',dpi=170);plt.close(fig)
    pd.DataFrame(rows).to_csv(OUT/'cycle_fft_amplitudes.csv',index=False,encoding='utf-8-sig')
    fig,axes=plt.subplots(3,4,figsize=(18,10),sharex=True)
    fullrows=[]
    for col,folder in enumerate(sorted(ROOT.glob('20261121_*_90s_multicosine_*'))):
        trial=folder.name;q=pd.read_csv(folder/f'{trial}_quest.csv');force=pd.read_csv(folder/f'{trial}_combalance.csv')
        elapsed=10+np.arange(8000)/100
        pc=np.interp(elapsed,q.stimulus_elapsed_s,q.pc_mapped_time_s)
        ys=[np.interp(pc,q.pc_mapped_time_s,q.ap_velocity_m_s),np.interp(pc,force.pc_estimated_sample_time_s,force.sensor_Fy_N),np.interp(pc,force.pc_estimated_sample_time_s,force.sensor_Mx_Nm)]
        for row,y in enumerate(ys):
            # One continuous 80s record; only the overall mean is removed.
            y=y-y.mean();n=len(y);f=np.fft.rfftfreq(n,.01);amp=2*abs(np.fft.rfft(y))/n
            amp[0]/=2;amp[-1]/=2;keep=(f>0)&(f<=3)
            ax=axes[row,col];ax.plot(f[keep],amp[keep],'.-',color='#222222',ms=2,lw=.8)
            for ff in np.arange(.05,1.500001,.05):ax.axvline(ff,color='#c43c39',alpha=.35,lw=.7)
            ax.set(yscale='log',xlim=(0,3));ax.grid(axis='y',alpha=.2)
            if col==0:ax.set_ylabel(labels[row])
            for ff,v in zip(f[keep],amp[keep]):fullrows.append(dict(trial=trial,channel=channels[row],frequency_hz=ff,one_sided_amplitude=v))
        axes[0,col].set_title(trial.replace('20261121_','').replace('_90s_multicosine_',' trial '));axes[-1,col].set_xlabel('Frequency [Hz]')
    fig.suptitle('Continuous 80s FFT (0.0125 Hz bins, rectangular window); red lines: stimulated frequencies')
    fig.tight_layout();fig.savefig(OUT/'all_trials_fft_80s.png',dpi=170);plt.close(fig)
    pd.DataFrame(fullrows).to_csv(OUT/'fft_80s_amplitudes.csv',index=False,encoding='utf-8-sig')

if __name__=='__main__':main()
