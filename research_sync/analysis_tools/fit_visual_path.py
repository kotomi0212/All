"""Conditional Tv/gamma identification: velocity -> closed-loop COM acceleration."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_periodic_response import write

G=9.80665
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def base_response(f,m,h,anthro,dynamics):
    mb=.978*m;j,l={'male':(.0473,.513),'female':(.0481,.504)}[anthro]
    lb=l*h;jb=j*m*h*h;A=jb/lb+mb*lb
    kp=10.53*A;kd=3.97*A;s=2j*np.pi*np.asarray(f);delay=np.exp(-.155*s)
    denominator=(kp+kd*s)*delay-mb*G
    if dynamics=='full':denominator+=A*s*s
    # Body output -> whole-body output, fixed feet. No platform inertial input.
    return (mb/m)*s*s*(-kp*delay)/denominator
def profile(f,y,base,mode):
    def conditional(t):
        b=base/(1+2j*np.pi*f*t)
        if mode=='gain_only':v=abs(b);gamma=max(float(v@abs(y)/(v@v)),0);err=gamma*v-abs(y)
        else:gamma=max(float(np.vdot(b,y).real/np.vdot(b,b).real),0);err=gamma*b-y
        return float(np.sum(abs(err)**2)),gamma
    grid=np.logspace(-6,6,900);values=np.array([conditional(t) for t in grid]);i=int(np.argmin(values[:,0]))
    lo=np.log(grid[max(i-1,0)]);hi=np.log(grid[min(i+1,len(grid)-1)])
    opt=minimize_scalar(lambda z:conditional(np.exp(z))[0],bounds=(lo,hi),method='bounded')
    tv=float(np.exp(opt.x));sse,gamma=conditional(tv);status='finite'
    zero_sse,zero_gamma=conditional(0)
    b=base/(2j*np.pi*f)
    if mode=='gain_only':rate=max(float(abs(b)@abs(y)/(abs(b)@abs(b))),0);infsse=float(np.sum((rate*abs(b)-abs(y))**2))
    else:rate=max(float(np.vdot(b,y).real/np.vdot(b,b).real),0);infsse=float(np.sum(abs(rate*b-y)**2))
    tolerance=1e-9*max(float(np.sum(abs(y)**2)),1e-30)
    if zero_sse<=sse+tolerance:tv=0.;gamma=zero_gamma;sse=zero_sse;status='T_zero_limit'
    if infsse<=sse+tolerance:tv=np.inf;gamma=np.inf;sse=infsse;status='T_infinite_limit_only_gamma_over_T_identified'
    if status=='finite':pred=gamma*base/(1+2j*np.pi*f*tv)
    elif tv==0:pred=gamma*base
    else:pred=rate*base/(2j*np.pi*f)
    close=grid[values[:,0]<=sse*1.05+tolerance]
    phase=float(np.sqrt(np.mean(np.angle(pred*np.conj(y))**2)))*180/np.pi if np.any(abs(pred)>0) else np.nan
    if not np.any(abs(pred)>0):
        tv=np.nan;gamma=0.;status='no_positive_gamma_fit_T_unidentifiable'
    gain_sse=float(np.sum((abs(pred)-abs(y))**2));sst=float(np.sum((abs(y)-np.mean(abs(y)))**2))
    summary=dict(Tv_s=tv,gamma_v_s=gamma,gamma_over_T_if_infinite=rate if np.isinf(tv) else '',status=status,SSE=sse,
                 gain_R2_percent=100*(1-gain_sse/sst) if sst>0 else np.nan,phase_RMSE_deg=phase,
                 T_5percent_SSE_low_s=float(close.min()) if len(close) else '',T_5percent_SSE_high_s=float(close.max()) if len(close) else '')
    return summary,grid,values,pred

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    analysis=args.analysis.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    manifest=read(analysis/'trial_manifest.csv');frf=read(analysis/'frequency_response_all.csv');cycles=read(analysis/'cycle_coefficients_all.csv')
    rows=[];profiles=[];points=[];loo=[]
    for m in manifest:
        tag=m['output'];folder=Path(m['source']);meta=json.loads((folder/(folder.name+'_metadata.json')).read_text(encoding='utf-8-sig'))
        mass=float(m['mass_kg']);height=meta['participant']['height_cm']/100
        rr=[r for r in frf if r['recording_tag']==tag and r['output']=='COM_acceleration_Fy_over_M']
        f=np.array([float(r['frequency_hz']) for r in rr]);y=np.array([complex(float(r['H_real']),float(r['H_imag'])) for r in rr]);mask=f<=.700000001
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        axes[0,0].loglog(f,abs(y),'ko',label='measured FRF');axes[1,0].semilogx(f,np.angle(y,deg=True),'ko',label='measured')
        for anthro in ['male','female']:
          for dynamics in ['full','reduced']:
           for mode in ['gain_only','complex']:
            for sign in ([1,-1] if mode=='complex' else [1]):
                base=base_response(f,mass,height,anthro,dynamics);obs=sign*y
                fit,grid,val,pred=profile(f[mask],obs[mask],base[mask],mode)
                ident=dict(participant=m['participant'],trial=m['trial_number'],condition=m['condition'],recording_tag=tag,
                           anthropometry_assumption=anthro,body_dynamics=dynamics,fit=mode,force_sign=sign,**fit)
                rows.append(ident)
                if anthro=='male' and dynamics=='full' and sign==1:
                    # Main figure: existing default male coefficients are conditional assumptions.
                    if np.isinf(fit['Tv_s']):densepred=fit['gamma_over_T_if_infinite']*base/(2j*np.pi*f)
                    else:densepred=fit['gamma_v_s']*base/(1+2j*np.pi*f*fit['Tv_s'])
                    dense=np.geomspace(f.min(),f.max(),800)
                    densebase=base_response(dense,mass,height,anthro,dynamics)
                    modelcurve=fit['gamma_over_T_if_infinite']*densebase/(2j*np.pi*dense) if np.isinf(fit['Tv_s']) else fit['gamma_v_s']*densebase/(1+2j*np.pi*dense*fit['Tv_s'])
                    color='#0072B2' if mode=='gain_only' else '#d62728'
                    axes[0,0].loglog(dense,abs(modelcurve),'-',color=color,label=mode)
                    axes[1,0].semilogx(dense,np.angle(modelcurve,deg=True),'-',color=color)
                    axes[0,1].semilogx(grid,val[:,0]/max(fit['SSE'],1e-30),color=color,label=mode)
                    axes[1,1].loglog(grid,np.maximum(val[:,1],1e-20),color=color,label=mode)
                    for t,(e,g) in zip(grid,val):profiles.append(dict(recording_tag=tag,fit=mode,Tv_s=t,gamma_v_s=g,SSE=e,SSE_over_best=e/max(fit['SSE'],1e-30)))
                    for hz,ob,pr in zip(f,obs,densepred):points.append(dict(recording_tag=tag,fit=mode,frequency_hz=hz,used_for_fit=hz<=.70000001,observed_gain=abs(ob),model_gain=abs(pr),observed_phase_deg=np.angle(ob,deg=True),model_phase_deg=np.angle(pr,deg=True)))
                    # Leave-one-cycle-out sensitivity, not confidence intervals.
                    cr=[r for r in cycles if r['recording_tag']==tag and r['output']=='COM_acceleration_Fy_over_M']
                    for excluded in range(1,5):
                        estimates=[]
                        for hz in f:
                            selected=[r for r in cr if int(r['cycle'])!=excluded and np.isclose(float(r['frequency_hz']),hz)]
                            U=np.array([complex(float(r['U_real']),float(r['U_imag'])) for r in selected]);Y=np.array([complex(float(r['Y_real']),float(r['Y_imag'])) for r in selected]);estimates.append(np.vdot(U,Y)/np.vdot(U,U).real)
                        fit_loo,*_=profile(f[mask],np.array(estimates)[mask],base[mask],mode)
                        loo.append(dict(recording_tag=tag,fit=mode,excluded_cycle=excluded,**fit_loo))
        axes[0,0].set(ylabel='COM acceleration / background velocity (1/s)');axes[1,0].set(xlabel='Frequency (Hz)',ylabel='Phase (deg)',ylim=(-190,190))
        for ax in axes[:,0]:ax.axvspan(.7,1.5,color='gray',alpha=.1);ax.grid(which='both',alpha=.2)
        axes[0,0].legend(fontsize=8);axes[0,1].set(xlabel='Tv (s)',ylabel='SSE / best SSE',ylim=(.98,4));axes[0,1].legend(fontsize=8)
        axes[1,1].set(xlabel='Tv (s)',ylabel='Re-estimated gamma (s)')
        for ax in axes[:,1]:ax.grid(which='both',alpha=.2)
        fig.suptitle(m['participant']+' trial '+m['trial'] if 'trial' in m else m['participant']+' T'+m['trial_number']+' | '+m['condition']+'\nConditional full closed-loop model; fixed controller, assumed body coefficients')
        fig.savefig(out/(tag+'_identification.png'),dpi=150);plt.close(fig)
    write(out/'all_conditional_fits.csv',rows);write(out/'Tv_profiles.csv',profiles);write(out/'model_vs_measured.csv',points);write(out/'leave_one_cycle_out.csv',loo)
    main=[r for r in rows if r['anthropometry_assumption']=='male' and r['body_dynamics']=='full' and r['force_sign']==1]
    write(out/'main_conditional_fits.csv',main)
    lines=['# Tv・gamma_vの条件付き同定','',
    '全5試行を別々に同定。入力=記録された背景AP速度、出力=Fy/M。0.05–0.70 Hzの7点を使用。0.9–1.5 Hzは図に表示のみ。時間シフトや追加遅延の最適化なし。',
    '主結果は身体慣性を含む閉ループ。A=Jb/lb+mb*lb、kp=10.53A、kd=3.97A、td=.155 sを既存コードから固定。男性係数Jb=.0473MH², lb=.513H, mb=.978Mを仮定し、女性係数での感度も計算。これらは個人別校正値ではない。',
    '式(8)は(Tv*s+1) X_target=-gamma_v V_background。固定床につき支持面加速度項0。Gx=-kp*gamma_v*exp(-td*s)/[(Tv*s+1){A*s²+(kp+kd*s)exp(-td*s)-mb*g}]。全身COM加速度/速度のモデルは(mb/M)*s²*Gx（足部静止）。gamma_vはこの定義では秒。一次遅れ単体をCOMに直接フィットしていない。',
    'gain_only: 線形ゲイン残差の最小二乗。complex: FRFの実部・虚部残差の最小二乗。各周波数は等重み。各Tvでgamma>=0を解析的に最適化。Tv=1e-6–1e6秒を走査しT=0とT→∞も比較。無限大極限ではgamma/Tのみ求まりTvとgammaの個別値を報告しない。',
    '', '|匿名被験者|条件|試行|損失|Tv s|gamma s|状態|ゲインR² %|位相RMSE °|','|---|---|---|---|---:|---:|---|---:|---:|']
    for r in main:lines.append(f"|{r['participant']}|{r['condition']}|{r['trial']}|{r['fit']}|{r['Tv_s']:.4g}|{r['gamma_v_s']:.4g}|{r['status']}|{r['gain_R2_percent']:.1f}|{r['phase_RMSE_deg']:.1f}|")
    lines+=['','SSEが最小値から5%以内のTv範囲は形状の感度指標であり信頼区間ではない。leave_one_cycle_out.csvは4周期のどれかを外した感度で、独立試行による検証でもない。',
    'all_conditional_fits.csvには、既存の身体慣性省略モデルとの比較、男女係数、複素フィットの力符号反転を含む。符号は良く合う方を正しいと選ばない。主結果は既存Fy/M規約。',
    '結果は仮定付きであり、原著の個人別Tv/gamma_vとして確定したものではない。低S/N・周期間ばらつき、未校正のAP符号・物理遅延が残る。以前の波形解析と併せて解釈する。元データは読み取りのみ。',
    '', '```powershell', '.\\.venv\\Scripts\\python.exe -X utf8 research_sync/analysis_tools/fit_visual_path.py --analysis research_sync/results/reanalysis_20261125_three_subjects --output research_sync/results/NEW_IDENTIFICATION', '```']
    (out/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8');print(json.dumps(main,indent=2))

if __name__=='__main__':main()
