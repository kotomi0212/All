"""Sonobe & Inoue 2019 Eq.(3), stationary plate; provisional axis mapping.

Run: python 261121/refit_sonobe2019.py
Preserves the previous reduced balance model and fitting band for comparison.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from fit_direct_com_acceleration import complex_amplitude, wrapped_phase_error
from compare_com_output_models import displacement_model
from fit_visual_acceleration_first_order import model as first_order_model, fit as first_order_fit
from fit_visual_acceleration_first_order_gain_only import fit_gain_only as first_order_gain_fit

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'sonobe2019_refit'
F = np.arange(.05, 1.500001, .05)
MASK = F <= .70000001
G = 9.80665

def fit(obs, pars, mode):
    def residual(z):
        pred = displacement_model(F[MASK], *np.exp(z), **pars)
        if mode == 'gain_only':
            return (abs(pred)-abs(obs[MASK])) / np.mean(abs(obs[MASK]))
        return np.r_[np.log(abs(pred)/abs(obs[MASK])), wrapped_phase_error(pred, obs[MASK])]
    results = [least_squares(residual, np.log([t,.03]), bounds=(np.log([.001,1e-8]),np.log([1000,100])), max_nfev=10000) for t in (.01,.3,3,100)]
    best = min(results, key=lambda r: np.sum(r.fun**2))
    tv, gamma = np.exp(best.x)
    pred = displacement_model(F[MASK], tv, gamma, **pars)
    r2 = 100*(1-np.sum((abs(pred)-abs(obs[MASK]))**2)/np.sum((abs(obs[MASK])-np.mean(abs(obs[MASK])))**2))
    return tv,gamma,r2,float(np.sqrt(np.mean((20*np.log10(abs(pred)/abs(obs[MASK])))**2))),float(np.sqrt(np.mean(np.degrees(wrapped_phase_error(pred,obs[MASK]))**2)))

def main():
    OUT.mkdir(exist_ok=True)
    rows, frfs, checks, first_rows, sensitivity = [], [], [], [], []
    fig, axes = plt.subplots(2,4,figsize=(16,7),sharex=True)
    for col,folder in enumerate(sorted(ROOT.glob('20261121_*_90s_multicosine_*'))):
        stem=folder.name
        q=pd.read_csv(folder/f'{stem}_quest.csv'); d=pd.read_csv(folder/f'{stem}_combalance.csv')
        H=float(q.height_cm.iloc[0])/100; M=float(q.body_mass_kg.iloc[0]); mb=.978*M
        qt=q.pc_mapped_time_s.to_numpy(); ft=d.pc_estimated_sample_time_s.to_numpy()
        start,end=np.interp([10,90],q.stimulus_elapsed_s,qt)
        keep=(ft>=start)&(ft<end); t=ft[keep]-start
        u=np.interp(ft[keep],qt,q.ap_velocity_m_s)
        rx=d.loc[keep,'sensor_Fy_N'].to_numpy(); ny=-d.loc[keep,'sensor_Mx_Nm'].to_numpy()
        ua=complex_amplitude(t,u,F)
        visual_acc=np.interp(ft[keep],qt,q.ap_acceleration_m_s2)
        visual_acc_amp=complex_amplitude(t,visual_acc,F)
        series=pd.DataFrame({'time_s':t,'vr_velocity_m_s':u,'Rx_assumed_N':rx,'Ny_assumed_Nm':ny})
        for sex,jratio,lratio in [('male',.0473,.513),('female',.0481,.504)]:
            jb=jratio*M*H**2; lb=lratio*H; lf=.039*H
            c=jb/(mb*lb)+lf+lb
            # Eq.(3): Mg*xp is represented by sagittal plate moment Ny.
            # Rx is force applied TO plate in the paper. Hardware sign remains unvalidated.
            x=(ny-c*rx)/(mb*G)
            pars={'body_mass':mb,'inertia_term':jb/lb+mb*lb}
            obs=complex_amplitude(t,x,F)/ua
            # Also reproduce the other previous model: visual acceleration -> COM acceleration.
            accel_frf=(1j*2*np.pi*F)**2*obs*ua/visual_acc_amp
            for method in ['gain_only','gain_phase']:
                if method=='gain_only':
                    tt,gg,rr,ee,_=first_order_gain_fit(accel_frf)
                    pp=first_order_model(F[MASK],tt,gg)
                    pe=float(np.sqrt(np.mean(np.degrees(wrapped_phase_error(pp,accel_frf[MASK]))**2)))
                else:
                    tt,gg,ee,pe=first_order_fit(accel_frf);rr=float('nan')
                first_rows.append(dict(trial=stem,anthropometry=sex,fit_method=method,T_v_s=tt,gamma_v=gg,inverse_gamma_v=1/gg,gain_R2_percent=rr,gain_rmse_db=ee,phase_rmse_deg=pe,T_v_at_bound=bool(tt<.00101 or tt>999)))
            # Sensor moment sign has conflicting assumptions in older code.
            # Record both signs, without choosing the one giving a better fit.
            for moment_sign in [-1,1]:
                xs=(moment_sign*d.loc[keep,'sensor_Mx_Nm'].to_numpy()-c*rx)/(mb*G)
                os=complex_amplitude(t,xs,F)/ua
                tt,gg,rr,ee,pe=fit(os,pars,'gain_only')
                sensitivity.append(dict(trial=stem,anthropometry=sex,moment_sign=moment_sign,T_v_s=tt,gamma_v_s=gg,gain_R2_percent=rr,phase_rmse_deg=pe))
            series[f'body_COM_{sex}_m']=x
            for mode in ['gain_only','gain_phase']:
                tv,gamma,r2,gr,pr=fit(obs,pars,mode)
                pred=displacement_model(F,tv,gamma,**pars)
                rows.append(dict(trial=stem,subject=q.subject.iloc[0],trial_number=int(q.trial_number.iloc[0]),anthropometry=sex,fit_method=mode,T_v_s=tv,gamma_v_s=gamma,inverse_gamma_v_per_s=1/gamma,W_v=gamma/tv,gain_R2_percent=r2,gain_rmse_db=gr,phase_rmse_deg=pr,T_v_at_bound=bool(tv<.00101 or tv>999),height_m=H,mass_kg=M,mass_status=d.body_mass_status.iloc[0],axis_status='provisional: Rx=Fy, Ny=-Mx',samples=len(t)))
                for f,o,p in zip(F,obs,pred):
                    frfs.append(dict(trial=stem,anthropometry=sex,fit_method=mode,frequency_hz=f,used_for_fit=bool(f<=.70000001),measured_gain_s=abs(o),model_gain_s=abs(p),measured_phase_deg=np.angle(o,deg=True),model_phase_deg=np.angle(p,deg=True)))
                if sex=='male':
                    dense=np.linspace(.05,1.5,500); pp=displacement_model(dense,tv,gamma,**pars)
                    axes[0,col].plot(dense,abs(pp),label=mode)
                    axes[1,col].plot(dense,np.angle(pp,deg=True))
            if sex=='male':
                axes[0,col].scatter(F,abs(obs),s=15,color='black',label='Eq.(3) data')
                axes[1,col].scatter(F,np.angle(obs,deg=True),s=15,color='black')
            # Newton acceleration from paper is -Rx/mb; compare, without enforcing equality.
            ratio=(1j*2*np.pi*F)**2*complex_amplitude(t,x,F)/complex_amplitude(t,-rx/mb,F)
            for f,r in zip(F,ratio):
                checks.append(dict(trial=stem,anthropometry=sex,frequency_hz=f,acceleration_ratio_gain=abs(r),acceleration_ratio_phase_deg=np.angle(r,deg=True)))
        series.to_csv(OUT/f'{stem}_com.csv',index=False,encoding='utf-8-sig')
        axes[0,col].set_title(f'{q.subject.iloc[0]} trial {q.trial_number.iloc[0]}')
        for ax in axes[:,col]:
            ax.axvspan(.7,1.5,color='gray',alpha=.12); ax.grid(alpha=.25); ax.set_xlim(.04,1.5)
        axes[0,col].set_yscale('log'); axes[1,col].set_ylim(-190,190); axes[1,col].set_xlabel('Frequency [Hz]')
    axes[0,0].set_ylabel('COM displacement / VR velocity [s]'); axes[1,0].set_ylabel('Phase [deg]')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Sonobe 2019 Eq.(3): male parameters; provisional force/moment signs')
    fig.tight_layout(); fig.savefig(OUT/'displacement_gain_phase.png',dpi=180); plt.close(fig)
    results=pd.DataFrame(rows); results.to_csv(OUT/'identified_parameters.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(frfs).to_csv(OUT/'frequency_response.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(checks).to_csv(OUT/'acceleration_consistency.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(first_rows).to_csv(OUT/'first_order_acceleration_parameters.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(sensitivity).to_csv(OUT/'moment_sign_sensitivity.csv',index=False,encoding='utf-8-sig')
    assert np.isfinite(results[['T_v_s','gamma_v_s','gain_rmse_db','phase_rmse_deg']]).all().all()
    # Dimensional/algebraic check of moment representation with synthetic COP.
    xp=.01; assert np.isclose((mb*G*xp)/(mb*G),xp)
    print(results[['subject','trial_number','anthropometry','fit_method','T_v_s','gamma_v_s','gain_R2_percent','phase_rmse_deg','T_v_at_bound']].to_string(index=False))

if __name__=='__main__':
    main()
