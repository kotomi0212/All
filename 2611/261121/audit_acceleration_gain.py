"""Audit measured gain construction without fitting a visual response model."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from refit_sonobe2019_acceleration import ROOT,F,complex_amplitude

OUT=ROOT/'acceleration_gain_audit'

def joint_amplitude(t,y):
    tc=t-np.mean(t)
    design=np.column_stack([np.ones(len(t)),tc]+[v for f in F for v in (np.cos(2*np.pi*f*t),np.sin(2*np.pi*f*t))])
    b=np.linalg.lstsq(design,y,rcond=None)[0]
    return b[2::2]-1j*b[3::2]

def main():
    OUT.mkdir(exist_ok=True);rows=[];summaries=[]
    fig,axes=plt.subplots(5,4,figsize=(17,15),sharex=True)
    fig2,axs=plt.subplots(3,4,figsize=(17,10),sharex=True)
    for col,folder in enumerate(sorted(ROOT.glob('20261121_*_90s_multicosine_*'))):
        stem=folder.name;q=pd.read_csv(folder/f'{stem}_quest.csv');d=pd.read_csv(folder/f'{stem}_combalance.csv')
        H=float(q.height_cm.iloc[0])/100;M=float(q.body_mass_kg.iloc[0]);mb=.978*M
        qt=q.pc_mapped_time_s.to_numpy();ft=d.pc_estimated_sample_time_s.to_numpy()
        start,end=np.interp([10,90],q.stimulus_elapsed_s,qt);keep=(ft>=start)&(ft<end);t=ft[keep]-start
        v=np.interp(ft[keep],qt,q.ap_velocity_m_s);u=np.interp(ft[keep],qt,q.ap_displacement_m)
        rx=d.loc[keep,'sensor_Fy_N'].to_numpy();ny=-d.loc[keep,'sensor_Mx_Nm'].to_numpy()
        jb=.0473*M*H**2;lb=.513*H;c=jb/(mb*lb)+.039*H+lb
        x=(ny-c*rx)/(mb*9.80665)
        va=complex_amplitude(t,v,F);xa=complex_amplitude(t,x,F);ua=complex_amplitude(t,u,F)
        ra=complex_amplitude(t,rx,F);na=complex_amplitude(t,ny,F)
        ratio=abs(xa/va);acc=(2*np.pi*F)**2*abs(xa);gain=acc/abs(va)
        vj=joint_amplitude(t,v);xj=joint_amplitude(t,x)
        # Four independent 20-second blocks; retain global stimulus phase.
        blocks=[]
        for k in range(4):
            m=(t>=20*k)&(t<20*(k+1))
            blocks.append(joint_amplitude(t[m],x[m])/joint_amplitude(t[m],v[m]))
        blocks=np.array(blocks)
        # Repeatability indicator only: not magnitude-squared coherence or a noise floor.
        repeat=abs(np.mean(blocks,axis=0))**2/np.mean(abs(blocks)**2,axis=0)
        fmask=F<=.7000001
        slopes={}
        for name,y in [('VR_velocity',abs(va)),('COM_displacement',abs(xa)),('displacement_gain',ratio),('acceleration_gain',gain)]:
            slopes[name]=np.polyfit(np.log10(F[fmask]),np.log10(y[fmask]),1)[0]
        summaries.append(dict(trial=stem,input_velocity_max_min_ratio=max(abs(va))/min(abs(va)),input_velocity_log_slope=slopes['VR_velocity'],COM_displacement_log_slope=slopes['COM_displacement'],displacement_gain_log_slope=slopes['displacement_gain'],acceleration_gain_log_slope=slopes['acceleration_gain'],median_repeatability_fitband=np.median(repeat[fmask]),median_joint_vs_old_gain_relative_difference=np.median(abs(abs(xj/vj)/ratio-1)[fmask])))
        for i,f in enumerate(F):
            rows.append(dict(trial=stem,frequency_hz=f,fit_band=bool(f<=.7000001),VR_displacement_amplitude_m=abs(ua[i]),VR_velocity_amplitude_m_s=abs(va[i]),COM_displacement_amplitude_m=abs(xa[i]),displacement_gain_s=ratio[i],acceleration_multiplier_per_s2=(2*np.pi*f)**2,COM_acceleration_amplitude_m_s2=acc[i],acceleration_gain_per_s=gain[i],force_contribution_to_COM_amplitude_m=abs(c*ra[i]/(mb*9.80665)),moment_contribution_to_COM_amplitude_m=abs(na[i]/(mb*9.80665)),joint_regression_acceleration_gain_per_s=(2*np.pi*f)**2*abs(xj[i]/vj[i]),block_repeatability=repeat[i],block_gain_amplitude_CV=np.std(abs(blocks[:,i]))/np.mean(abs(blocks[:,i])),**{f'block{k+1}_acceleration_gain_per_s':abs(blocks[k,i])*(2*np.pi*f)**2 for k in range(4)}))
        yy=[abs(va),abs(xa)*1000,ratio,acc,gain]
        for row,y in enumerate(yy):
            axes[row,col].plot(F,y,'o-',ms=3,lw=.8,color='#222222')
            axes[row,col].set_yscale('log');axes[row,col].axvspan(.7,1.5,color='gray',alpha=.12);axes[row,col].grid(alpha=.25)
        axs[0,col].plot(F,abs(na/(mb*9.80665))*1000,'o-',ms=3,label='Moment contribution')
        axs[0,col].plot(F,abs(c*ra/(mb*9.80665))*1000,'o-',ms=3,label='Force contribution')
        axs[0,col].plot(F,abs(xa)*1000,'k.-',label='Combined COM')
        axs[0,col].set_yscale('log')
        axs[1,col].plot(F,gain,'k.-',label='Previous extraction')
        axs[1,col].plot(F,(2*np.pi*F)**2*abs(xj/vj),'x--',label='Joint extraction')
        axs[1,col].set_yscale('log')
        for k in range(4):axs[2,col].plot(F,abs(blocks[k])*(2*np.pi*F)**2,'.-',alpha=.65,lw=.8,label=f'20s block {k+1}')
        axs[2,col].set_yscale('log')
        for ax in axs[:,col]:ax.grid(alpha=.25);ax.axvspan(.7,1.5,color='gray',alpha=.12)
        title=f'{q.subject.iloc[0]} trial {int(q.trial_number.iloc[0])}'
        axes[0,col].set_title(title);axs[0,col].set_title(title)
        axes[-1,col].set_xlabel('Frequency [Hz]');axs[-1,col].set_xlabel('Frequency [Hz]')
        assert np.allclose(gain,ratio*(2*np.pi*F)**2)
    labels=['1. VR velocity amplitude [m/s]','2. COM displacement amplitude [mm]','3. COM displacement / VR velocity [s]','4. COM acceleration amplitude [m/s2]','5. COM acceleration / VR velocity [1/s]']
    for ax,label in zip(axes[:,0],labels):ax.set_ylabel(label)
    for ax,label in zip(axs[:,0],['COM contributions [mm]','Acceleration gain [1/s]','20s block acceleration gain [1/s]']):ax.set_ylabel(label);ax.legend(fontsize=7)
    fig.suptitle('Data construction audit: male Sonobe parameters; no visual response model fitted')
    fig.tight_layout();fig.savefig(OUT/'gain_construction_stages.png',dpi=160)
    fig2.suptitle('Force/moment contributions, extraction sensitivity, and four-cycle repeatability')
    fig2.tight_layout();fig2.savefig(OUT/'processing_and_repeatability.png',dpi=160)
    pd.DataFrame(rows).to_csv(OUT/'frequency_audit.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(summaries).to_csv(OUT/'summary.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(summaries).to_string(index=False))
    print(pd.DataFrame(rows).query('frequency_hz < 0.101 or (frequency_hz > 0.499 and frequency_hz < 0.501)').to_string(index=False))

if __name__=='__main__':main()
