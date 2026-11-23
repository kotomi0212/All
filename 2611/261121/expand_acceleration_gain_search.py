"""Expanded gain-only search with exact conditional optimization of gamma."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from refit_sonobe2019_acceleration import ROOT, F, MASK, model

OUT=ROOT/'sonobe2019_acceleration_wide_search'
T_MIN,T_MAX=1e-8,1e6
GAMMA_MIN,GAMMA_MAX=1e-12,1e6

def main():
    OUT.mkdir(exist_ok=True)
    source=ROOT/'sonobe2019_acceleration_refit'
    data=pd.read_csv(source/'frequency_response.csv')
    old=pd.read_csv(source/'identified_parameters.csv')
    summaries=[];profiles=[];responses=[]
    fig,axes=plt.subplots(2,4,figsize=(16,7))
    for col,trial in enumerate(data.trial.unique()):
        q=pd.read_csv(ROOT/trial/f'{trial}_quest.csv',nrows=1)
        H=float(q.height_cm.iloc[0])/100;M=float(q.body_mass_kg.iloc[0]);mb=.978*M
        for sex,jr,lr in [('male',.0473,.513),('female',.0481,.504)]:
            d=data[(data.trial==trial)&(data.anthropometry==sex)&(data.fit_method=='gain_only')].sort_values('frequency_hz')
            f=d.frequency_hz.to_numpy();y=d.measured_gain_per_s.to_numpy();mask=f<=.70000001
            lb=lr*H;jb=jr*M*H**2;pars={'body_mass':mb,'inertia_term':jb/lb+mb*lb}
            def conditional(t):
                b=abs(model(f[mask],t,1.,pars))
                gamma=np.clip(np.dot(b,y[mask])/np.dot(b,b),GAMMA_MIN,GAMMA_MAX)
                return float(np.sum((gamma*b-y[mask])**2)),float(gamma)
            grid=np.linspace(np.log(T_MIN),np.log(T_MAX),2801)
            errors=np.array([conditional(np.exp(z))[0] for z in grid])
            k=int(np.argmin(errors))
            candidates=[(errors[0],T_MIN),(errors[-1],T_MAX)]
            # Refine all local minima after the global grid scan.
            for i in range(1,len(grid)-1):
                if errors[i]<errors[i-1] and errors[i]<=errors[i+1]:
                    opt=minimize_scalar(lambda z:conditional(np.exp(z))[0],bounds=(grid[i-1],grid[i+1]),method='bounded')
                    candidates.append((opt.fun,float(np.exp(opt.x))))
            sse,tv=min(candidates,key=lambda v:v[0])
            # Roundoff on the flat T->0 plateau must not imply an interior minimum.
            if errors[0]<=sse+1e-12*max(1.,sse):
                tv=T_MIN;sse=errors[0]
            _,gamma=conditional(tv)
            zero_sse,zero_gamma=conditional(0.)
            o=old[(old.trial==trial)&(old.anthropometry==sex)&(old.fit_method=='gain_only')].iloc[0]
            old_sse=conditional(float(o.T_v_s))[0]
            r2=100*(1-sse/np.sum((y[mask]-np.mean(y[mask]))**2))
            summaries.append(dict(trial=trial,anthropometry=sex,T_v_s=tv,gamma_v_s=gamma,gain_R2_percent=r2,T_at_lower_bound=bool(tv<=T_MIN*1.001),gamma_at_bound=bool(gamma<=GAMMA_MIN*1.001 or gamma>=GAMMA_MAX/.999),old_T_v_s=o.T_v_s,old_gamma_v_s=o.gamma_v_s,old_R2_percent=o.gain_R2_percent,SSE=sse,SSE_improvement_percent=100*(old_sse-sse)/old_sse,T_zero_limit_SSE=zero_sse,T_zero_limit_gamma=zero_gamma,T_min_s=T_MIN,T_max_s=T_MAX,gamma_min_s=GAMMA_MIN,gamma_max_s=GAMMA_MAX))
            assert sse<=old_sse*(1+1e-9)
            for z,e in zip(grid,errors):
                t=np.exp(z);profiles.append(dict(trial=trial,anthropometry=sex,T_v_s=t,gamma_opt_s=conditional(t)[1],SSE=e,SSE_ratio_to_best=e/sse))
            for ff,yy,pp in zip(f,y,abs(model(f,tv,gamma,pars))):
                responses.append(dict(trial=trial,anthropometry=sex,frequency_hz=ff,used_for_fit=bool(ff<=.70000001),measured_gain_per_s=yy,model_gain_per_s=pp))
            if sex=='male':
                dense=np.linspace(.05,1.5,600)
                axes[0,col].plot(dense,abs(model(dense,float(o.T_v_s),float(o.gamma_v_s),pars)),label='Previous bounds',lw=2,ls='--')
                axes[0,col].plot(dense,abs(model(dense,tv,gamma,pars)),label='Expanded bounds',lw=1.5)
                axes[0,col].scatter(f[mask],y[mask],s=17,c='black',label='Fit data')
                axes[0,col].scatter(f[~mask],y[~mask],s=17,facecolors='none',edgecolors='gray')
                axes[0,col].axvspan(.7,1.5,color='gray',alpha=.12)
                axes[0,col].set(yscale='log',xlabel='Frequency [Hz]',title=f'{q.subject.iloc[0]} trial {q.trial_number.iloc[0]}')
                axes[0,col].text(.04,.96,f'T = {tv:.1e} s (lower bound)\ngamma = {gamma:.5f} s\nR2 = {r2:.2f}%',transform=axes[0,col].transAxes,va='top',fontsize=9)
                axes[1,col].semilogx(np.exp(grid),errors/sse)
                axes[1,col].axvline(.001,c='gray',ls='--',label='Old T lower bound')
                axes[1,col].set(xlabel='T [s]; gamma optimized at each T',ylabel='SSE / minimum SSE',ylim=(.98,min(5,max(errors/sse)*1.02)))
        for ax in axes[:,col]:ax.grid(alpha=.25)
    axes[0,0].set_ylabel('COM acceleration / VR velocity [1/s]');axes[0,0].legend(fontsize=8,loc='lower right')
    fig.suptitle('Expanded gain-only search: T=1e-8..1e6 s, gamma=1e-12..1e6 s; male parameters')
    fig.tight_layout();fig.savefig(OUT/'gain_and_parameter_profile.png',dpi=180);plt.close(fig)
    pd.DataFrame(summaries).to_csv(OUT/'identified_parameters.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(profiles).to_csv(OUT/'T_profile.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(responses).to_csv(OUT/'frequency_response.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(summaries).to_string(index=False))

if __name__=='__main__':main()
