"""Quiet-standing comparison and noise-weighted conditional Tv/gamma fit."""
import argparse,csv,io,json,hashlib
from pathlib import Path
import numpy as np
from scipy.signal import welch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_periodic_response import harmonic,write
from fit_visual_path import base_response,profile,read

def main():
    p=argparse.ArgumentParser();p.add_argument('--static',nargs='+',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    frf=read(a.analysis/'frequency_response_all.csv');rep=read(a.analysis/'cycle_coefficients_all.csv');manifest=read(a.analysis/'trial_manifest.csv')
    F=np.array([.05,.1,.2,.3,.45,.6,.7,.9,1.2,1.5]);static_coeff=[];audit=[];bands=[]
    fig,ax=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for i,path in enumerate(a.static):
        raw=path.read_bytes();rows=list(csv.reader(io.StringIO(raw.decode('cp932'))));indices=[i for i,r in enumerate(rows) if 'rawFx' in r];h=rows[indices[-1]];meta={r[0]:r[1] for r in rows[:indices[0]] if len(r)==2}
        d=np.array([[float(v) for v in r] for r in rows[indices[-1]+1:] if r]);t=d[:,0]
        mass=float(meta['体重(kg)']);acc=d[:,h.index('Fx')]/mass;x=d[:,h.index('xb')]/1000
        if not np.allclose(np.diff(t),.01):raise ValueError('Static sampling')
        keep=(t>=5)&(t<55)
        for j,y in enumerate([acc,x]):
            f,psd=welch(y[keep],fs=100,nperseg=2000,noverlap=1000,detrend='linear',scaling='density')
            mask=(f>=.05)&(f<=2);cum=np.cumsum(psd[mask]);ff=f[mask]
            ax[j].loglog(f[mask],psd[mask],label=f'Static {i+1}')
            for q in [.5,.9,.95]:bands.append(dict(reference='Static '+str(i+1),output='Fx_over_M_acceleration' if j==0 else 'vendor_xb_position',analyzed_band_hz='0.05-2',power_fraction=q,frequency_at_cumulative_fraction_hz=float(ff[np.searchsorted(cum,q*cum[-1])]),note='descriptive, depends on output and band; not an exclusion cutoff'))
        for start in [5.,25.]:
            k=(t>=start)&(t<start+20);z=harmonic(t[k],acc[k],F);static_coeff.append(z)
        audit.append(dict(reference='Static '+str(i+1),path=str(path),sha256=hashlib.sha256(raw).hexdigest(),date=meta.get('試験開始日時'),condition=meta.get('実験の種類'),eyes=meta.get('開眼/閉眼'),mass_kg=mass,AP_column='Fx (OEM AP convention); lab axes vs VR unvalidated',position='vendor xb, processed; used as reference only',windows='5-25,25-45 s',static_trial_count=len(a.static),comparable='historical same-subject reference; unmatched day/equipment/processing'))
    for j,title in enumerate(['Static force-based COM acceleration','Static vendor COM position']):
        ax[j].set(title=title,xlim=(.05,2),xlabel='Frequency (Hz)',ylabel='PSD ((m/s²)²/Hz)' if j==0 else 'PSD (m²/Hz)');ax[j].legend();ax[j].grid(which='both',alpha=.2)
    fig.savefig(out/'static_spectra.png',dpi=150);plt.close(fig)
    static_coeff=np.array(static_coeff);baseline_power=np.mean(abs(static_coeff)**2,axis=0)
    comparisons=[];fits=[];prof=[]
    fig,axes=plt.subplots(3,1,figsize=(11,9),layout='constrained');seg_index=0
    for m in manifest:
        rr=[r for r in frf if r['recording_tag']==m['output'] and r['output']=='COM_acceleration_Fy_over_M'];R=[r for r in rep if r['recording_tag']==m['output'] and r['output']=='COM_acceleration_Fy_over_M']
        H=np.array([complex(float(r['H_real']),float(r['H_imag'])) for r in rr]);U=np.array([[complex(float(r['U_real']),float(r['U_imag'])) for r in R if int(r['cycle'])==k] for k in range(1,5)]);Y=np.array([[complex(float(r['Y_real']),float(r['Y_imag'])) for r in R if int(r['cycle'])==k] for k in range(1,5)])
        Z=Y/U;variance=np.var(Z,axis=0,ddof=1)/4
        config=json.loads((Path(m['source'])/(Path(m['source']).name+'_metadata.json')).read_text(encoding='utf-8-sig'));mass=float(m['mass_kg']);height=config['participant']['height_cm']/100
        base=base_response(F,mass,height,'male','full');mask=F<=.70000001
        weight_modes={'equal_weights':np.ones(10),'cycle_variance_weights':1/variance}
        # Static participant must match the named acquisition identity; never transfer to other participants.
        static_subjects=[path.name.split('_gridstatics_')[1].split('000')[0] for path in a.static]
        matched=config['participant']['subject'] in static_subjects
        if matched:weight_modes['historical_static_weights']=4*np.mean(abs(U)**2,axis=0)/baseline_power
        for mode,weights in weight_modes.items():
            scale=np.sqrt(weights[mask]/np.median(weights[mask]))
            fit,grid,values,pred=profile(F[mask],H[mask]*scale,base[mask]*scale,'complex')
            fits.append(dict(recording_tag=m['output'],participant=m['participant'],trial=m['trial_number'],condition=m['condition'],weighting=mode,**fit,note='R2/phase metrics in weighted space; compare prediction separately; no subtraction of static waveforms'))
            for tv,(sse,gamma) in zip(grid,values):prof.append(dict(recording_tag=m['output'],weighting=mode,Tv_s=tv,gamma_v_s=gamma,weighted_SSE=sse))
        if matched:
            k=seg_index;seg_index+=1
            coherent=abs(np.mean(Y,axis=0))**2
            finite_bias=np.var(Y,axis=0,ddof=1)/4
            corrected=np.maximum(coherent-finite_bias,0)
            axes[k].loglog(F,np.sqrt(coherent),'o-',label='Stimulus-aligned mean amplitude')
            axes[k].loglog(F,np.sqrt(baseline_power),'s--',label='Historical static RMS tone amplitude (20 s)')
            axes[k].loglog(F,np.where(corrected>0,np.sqrt(corrected),np.nan),'x:',label='Finite-cycle noise-bias adjusted mean amplitude')
            axes[k].set(title=m['participant']+' T'+m['trial_number']+' | '+m['condition'],xlabel='Frequency (Hz)',ylabel='COM acceleration (m/s²)');axes[k].grid(which='both',alpha=.2);axes[k].legend(fontsize=7)
            for i,f in enumerate(F):comparisons.append(dict(recording_tag=m['output'],frequency_hz=f,stim_aligned_mean_amplitude=np.sqrt(coherent[i]),static_rms_tone_amplitude=np.sqrt(baseline_power[i]),stim_vs_static_power_db=10*np.log10(coherent[i]/baseline_power[i]),bias_adjusted_coherent_power=corrected[i],interpretation='reference ratio, not validated SNR: baseline unmatched'))
    fig.suptitle('Historical static comparison | same 20 s harmonic estimator | not waveform subtraction')
    fig.savefig(out/'stimulus_vs_static.png',dpi=150);plt.close(fig)
    write(out/'static_inventory.csv',audit);write(out/'static_frequency_bands.csv',bands);write(out/'stimulus_static_comparison.csv',comparisons);write(out/'noise_weighted_fits.csv',fits);write(out/'weighted_Tv_profiles.csv',prof)
    lines=['# 自然動揺の影響評価と重み付き同定','',
    'static2試行（60秒・開眼安静立位、記録日時は2026/11/02）を確認。VR試行は11/24で日付・計測方向・装着条件が一致しないため、自然動揺の同時基準としては使えない。20秒×2区間/試行、計4区間を主試行と同じ多周波数回帰で比較。',
    '自然動揺は単一周波数のノイズではない。出力が位置か加速度かによってPSDの形も変わる。static_frequency_bands.csvの50/90/95%点は指定帯域内の記述値であり、一般的なカットオフではない。',
    '静止時の波形・振幅を刺激試行から直接引いていない。独立した揺れは位相が異なり、同じ周波数に刺激応答も含まれる。パワー差をとっても時間波形や刺激同期位相は復元できない。',
    'cycle_variance_weights: 各周期のFRF=Y/Uの複素分散/4を推定し、その逆数で実部・虚部の残差を重み付け。ばらつく周波数の寄与を下げるが、応答を新しく作り出す処理ではない。4周期で推定した分散自体にも不確実性がある。',
    'historical_static_weights: 同じ被験者のstaticだけを使用し、FRF分散≈静止時の20秒複素係数パワー/(4|U|²)とする探索的重み付け。別日の参考値であり主結果の補正とは断定しない。他の被験者へ流用しない。',
    'Tv/gammaの身体・制御係数は前回と同じ仮定（full、男性人口係数、kp=10.53A,kd=3.97A,td=.155s）。自然動揺への対処だけではこの前提や符号・同期・低S/Nの問題は解決しない。gain_R2は重み付け空間での値なので方法間の優劣比較に使わない。',
    '', '|記録|重み|Tv s|gamma s|状態|','|---|---|---:|---:|---|']
    for r in fits:lines.append(f"|{r['recording_tag']}|{r['weighting']}|{r['Tv_s']:.4g}|{r['gamma_v_s']:.4g}|{r['status']}|")
    lines+=['','参考文献： https://pmc.ncbi.nlm.nih.gov/articles/PMC4166606/ 、 https://pmc.ncbi.nlm.nih.gov/articles/PMC7237531/ 。刺激周波数成分の統計的検出・刺激との位相同期を評価する研究であり、特定帯域を自然動揺として削除する根拠にはしていない。']
    (out/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(static_trials=len(a.static),static_windows=len(static_coeff),fit_count=len(fits),matched_comparison_rows=len(comparisons)),indent=2))

if __name__=='__main__':main()
