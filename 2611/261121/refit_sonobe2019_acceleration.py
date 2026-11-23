"""Fit VR velocity -> acceleration derived from Sonobe 2019 COM displacement."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from refit_sonobe2019 import ROOT, F, MASK, complex_amplitude, wrapped_phase_error
from compare_com_output_models import displacement_model

OUT = ROOT / 'sonobe2019_acceleration_refit'

def model(f, tv, gamma, pars):
    return (1j*2*np.pi*f)**2 * displacement_model(f, tv, gamma, **pars)

def fit(obs, pars, mode):
    def residual(z):
        pred=model(F[MASK],*np.exp(z),pars)
        if mode=='gain_only':
            return (abs(pred)-abs(obs[MASK]))/np.mean(abs(obs[MASK]))
        return np.r_[np.log(abs(pred)/abs(obs[MASK])),wrapped_phase_error(pred,obs[MASK])]
    fits=[least_squares(residual,np.log([t,.03]),bounds=(np.log([.001,1e-8]),np.log([1000,100])),max_nfev=10000) for t in [.01,.3,3,100]]
    result=min(fits,key=lambda r:sum(r.fun**2))
    tv,gamma=np.exp(result.x);pred=model(F[MASK],tv,gamma,pars)
    r2=100*(1-sum((abs(pred)-abs(obs[MASK]))**2)/sum((abs(obs[MASK])-np.mean(abs(obs[MASK])))**2))
    gr=np.sqrt(np.mean((20*np.log10(abs(pred)/abs(obs[MASK])))**2))
    pr=np.sqrt(np.mean(np.degrees(wrapped_phase_error(pred,obs[MASK]))**2))
    return tv,gamma,r2,gr,pr

def main():
    OUT.mkdir(exist_ok=True);rows=[];values=[];harmonics=[]
    fig,axes=plt.subplots(2,4,figsize=(16,7),sharex=True)
    for col,folder in enumerate(sorted(ROOT.glob('20261121_*_90s_multicosine_*'))):
        stem=folder.name;q=pd.read_csv(folder/f'{stem}_quest.csv');d=pd.read_csv(folder/f'{stem}_combalance.csv')
        H=float(q.height_cm.iloc[0])/100;M=float(q.body_mass_kg.iloc[0]);mb=.978*M
        qt=q.pc_mapped_time_s.to_numpy();ft=d.pc_estimated_sample_time_s.to_numpy()
        start,end=np.interp([10,90],q.stimulus_elapsed_s,qt);keep=(ft>=start)&(ft<end);t=ft[keep]-start
        u=np.interp(ft[keep],qt,q.ap_velocity_m_s);ua=complex_amplitude(t,u,F)
        rx=d.loc[keep,'sensor_Fy_N'].to_numpy();ny=-d.loc[keep,'sensor_Mx_Nm'].to_numpy()
        for sex,jr,lr in [('male',.0473,.513),('female',.0481,.504)]:
            jb=jr*M*H**2;lb=lr*H;c=jb/(mb*lb)+.039*H+lb
            x=(ny-c*rx)/(mb*9.80665);xa=complex_amplitude(t,x,F);aa=(1j*2*np.pi*F)**2*xa;obs=aa/ua
            pars={'body_mass':mb,'inertia_term':jb/lb+mb*lb}
            for f,xx,a in zip(F,xa,aa):
                harmonics.append(dict(trial=stem,anthropometry=sex,frequency_hz=f,COM_displacement_real_m=xx.real,COM_displacement_imag_m=xx.imag,COM_acceleration_real_m_s2=a.real,COM_acceleration_imag_m_s2=a.imag))
            for mode in ['gain_only','gain_phase']:
                tv,gamma,r2,gr,pr=fit(obs,pars,mode);pred=model(F,tv,gamma,pars)
                rows.append(dict(trial=stem,subject=q.subject.iloc[0],trial_number=int(q.trial_number.iloc[0]),anthropometry=sex,fit_method=mode,input='VR_velocity_m_s',output='Sonobe_COM_acceleration_m_s2',T_v_s=tv,gamma_v_s=gamma,inverse_gamma_v_per_s=1/gamma,gain_R2_percent=r2,gain_rmse_db=gr,phase_rmse_deg=pr,T_v_at_bound=bool(tv<.00101 or tv>999),axis_status='provisional: Rx=Fy; Ny=-Mx'))
                for f,o,p in zip(F,obs,pred):
                    values.append(dict(trial=stem,anthropometry=sex,fit_method=mode,frequency_hz=f,used_for_fit=bool(f<=.70000001),measured_gain_per_s=abs(o),model_gain_per_s=abs(p),measured_phase_deg=np.angle(o,deg=True),model_phase_deg=np.angle(p,deg=True)))
                if sex=='male':
                    dense=np.linspace(.05,1.5,500);pp=model(dense,tv,gamma,pars)
                    axes[0,col].plot(dense,abs(pp),label=mode);axes[1,col].plot(dense,np.angle(pp,deg=True))
                if mode=='gain_phase':
                    # Log gain and phase residuals must be invariant to consistent s^2 output conversion.
                    px=displacement_model(F[MASK],tv,gamma,**pars);ox=xa[MASK]/ua[MASK]
                    assert np.allclose(np.log(abs(px)/abs(ox)),np.log(abs(pred[MASK])/abs(obs[MASK])))
                    assert np.allclose(wrapped_phase_error(px,ox),wrapped_phase_error(pred[MASK],obs[MASK]))
            if sex=='male':
                axes[0,col].scatter(F,abs(obs),s=15,c='black',label='COM acceleration data')
                axes[1,col].scatter(F,np.angle(obs,deg=True),s=15,c='black')
        axes[0,col].set_title(f'{q.subject.iloc[0]} trial {q.trial_number.iloc[0]}')
        for ax in axes[:,col]:
            ax.axvspan(.7,1.5,color='gray',alpha=.12);ax.grid(alpha=.25);ax.set_xlim(.04,1.5)
        axes[0,col].set_yscale('log');axes[1,col].set_ylim(-190,190);axes[1,col].set_xlabel('Frequency [Hz]')
    axes[0,0].set_ylabel('COM acceleration / VR velocity [1/s]');axes[1,0].set_ylabel('Phase [deg]');axes[0,0].legend(fontsize=8)
    fig.suptitle('Sonobe COM acceleration output: male parameters; provisional signs; fit 0.05-0.70 Hz')
    fig.tight_layout();fig.savefig(OUT/'acceleration_gain_phase.png',dpi=180);plt.close(fig)
    result=pd.DataFrame(rows);assert np.isfinite(result[['T_v_s','gamma_v_s','gain_R2_percent']]).all().all()
    result.to_csv(OUT/'identified_parameters.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(values).to_csv(OUT/'frequency_response.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(harmonics).to_csv(OUT/'COM_harmonic_acceleration.csv',index=False,encoding='utf-8-sig')
    print(result.to_string(index=False))

if __name__=='__main__':main()
