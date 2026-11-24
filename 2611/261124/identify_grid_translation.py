"""Unsynchronized 0.2/0.5 Hz translation: available HEAD output only.

This does not identify the COM/balance-controller parameters. Original force
plate CSVs are unavailable; exported head/modal acceleration must not replace COM.
Run: python 261124/identify_grid_translation.py
"""
from pathlib import Path
import hashlib
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

BASE = Path('C:/Users/atvma/gitseguchi/All/2610/sonobe_caluculation/outputs')
OUT = Path(__file__).resolve().parent / 'grid_translation_identification'
A = 0.0246  # displacement PEAK amplitude in m, supplied by user for both conditions
WINDOWS = [(5,35), (5,25), (10,30), (15,35), (10,50), (0,60)]


def harmonic(t, y, f):
    tc = t-t.mean()
    design = np.column_stack([np.ones(len(t)), tc, np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t)])
    coeff = np.linalg.lstsq(design, y, rcond=None)[0]
    trend = design[:,:2] @ np.linalg.lstsq(design[:,:2], y, rcond=None)[0]
    residual = y-design@coeff
    return float(np.hypot(*coeff[2:])), float(1-np.sum(residual**2)/np.sum((y-trend)**2))


def fit_gains(f, y):
    """Linear gain SSE, gamma optimized exactly; include T=0 and T=infinity."""
    w = 2*np.pi*np.asarray(f)
    y = np.asarray(y)
    def conditional(t):
        basis = 1/np.sqrt(1+(w*t)**2)
        gamma = np.dot(basis,y)/np.dot(basis,basis)
        return float(np.sum((gamma*basis-y)**2)), float(gamma)
    opt = minimize_scalar(lambda z: conditional(np.exp(z))[0], bounds=(-18,18), method='bounded')
    basis_inf = 1/w
    winf = np.dot(basis_inf,y)/np.dot(basis_inf,basis_inf)
    sinf = np.sum((winf*basis_inf-y)**2)
    choices = [(conditional(0)[0], 'T_zero', 0., conditional(0)[1], np.nan),
               (conditional(np.exp(opt.x))[0], 'finite', np.exp(opt.x), conditional(np.exp(opt.x))[1], np.nan),
               (sinf, 'T_infinite_limit', np.inf, np.inf, winf)]
    best = min(choices, key=lambda row: row[0])
    # Floating-point equality at large T must not be reported as finite identification.
    if sinf <= best[0]+1e-12*max(float(np.dot(y,y)),1e-30):
        best = choices[2]
    return dict(SSE=best[0], status=best[1], T_s=best[2], gamma_dimensionless=best[3], gamma_over_T_per_s=best[4]), conditional


def main():
    OUT.mkdir(exist_ok=True)
    # Check the actual estimator with offset, drift and arbitrary unknown phase.
    t = np.arange(3000)/100
    amplitude, r2 = harmonic(t, .1+.002*t+.007*np.cos(2*np.pi*.2*t+.83), .2)
    assert abs(amplitude-.007)<1e-12 and abs(r2-1)<1e-12
    synthetic = 0.04/np.sqrt(1+(2*np.pi*np.array([.2,.5])*1.2)**2)
    fit,_ = fit_gains([.2,.5],synthetic)
    assert abs(fit['T_s']-1.2)<1e-4 and abs(fit['gamma_dimensionless']-.04)<1e-6

    rows, provenance = [], []
    for p in sorted((BASE/'head_imu_20261026').glob('*gridtrans*waveform.csv')):
        d = pd.read_csv(p)
        f = .2 if 'trans02' in p.name else .5
        trial = p.stem.removesuffix('_waveform')
        joint = BASE/'joint_acceleration_20261101'/f'{trial}_joint_acceleration.csv'
        j = pd.read_csv(joint)
        assert np.allclose(j.time_s,d.time_s) and len(d)==6000
        # Modal sum is head acceleration, NOT COM acceleration.
        assert np.max(abs(j.AP_ankle_xix1_mm_s2+j.AP_hip_xix2_mm_s2-d.AP_axhd_mm_s2))<1e-6
        for source in [p,joint]:
            provenance.append(dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
        signals = {'head':d.AP_axhd_mm_s2, 'ankle_mode':j.AP_ankle_xix1_mm_s2, 'hip_mode':j.AP_hip_xix2_mm_s2}
        for start,end in WINDOWS:
            keep=(d.time_s>=start)&(d.time_s<end)
            t=d.time_s[keep].to_numpy()
            for name,values in signals.items():
                y=values[keep].to_numpy()/1000
                amp,explained=harmonic(t,y,f)
                # Descriptive peak only, not a replacement for the known input frequency.
                ff=np.fft.rfftfreq(len(y),.01); spec=abs(np.fft.rfft((y-y.mean())*np.hanning(len(y))))
                band=(ff>=.1)&(ff<=1.)
                peak=ff[band][np.argmax(spec[band])]
                rows.append(dict(trial=trial,output=name,frequency_Hz=f,start_s=start,end_s=end,
                    samples=len(t),input_displacement_peak_m=A,input_velocity_peak_m_s=A*2*np.pi*f,
                    input_acceleration_peak_m_s2=A*(2*np.pi*f)**2,
                    output_acceleration_peak_m_s2=amp,gain_acceleration_ratio=amp/(A*(2*np.pi*f)**2),
                    harmonic_partial_R2=explained,dominant_frequency_0p1_1_Hz=peak))
    data=pd.DataFrame(rows)
    assert len(data)==4*len(WINDOWS)*3 and np.isfinite(data.select_dtypes('number')).all().all()
    data.to_csv(OUT/'harmonic_amplitudes.csv',index=False)
    pd.DataFrame(provenance).to_csv(OUT/'provenance.csv',index=False)
    fits=[]; profiles=[]
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    primary=data[(data.output=='head')&(data.start_s==5)&(data.end_s==35)]
    for start,end in WINDOWS:
        block=data[(data.output=='head')&(data.start_s==start)&(data.end_s==end)]
        fit,conditional=fit_gains(block.frequency_Hz,block.gain_acceleration_ratio)
        means=block.groupby('frequency_Hz').gain_acceleration_ratio.mean()
        ratio=means.loc[.5]/means.loc[.2]
        fits.append(dict(output='head',start_s=start,end_s=end,selection='all_four_trials',gain_0p5_over_0p2=ratio,**fit))
        ts=np.logspace(-3,5,500)
        ss=np.array([conditional(t)[0] for t in ts])
        for tv,sse in zip(ts,ss):
            profiles.append(dict(start_s=start,end_s=end,T_s=tv,SSE=sse,SSE_over_minimum=sse/fit['SSE']))
        axes[1].semilogx(ts,ss/fit['SSE'],label=f'{start}-{end} s')
    for (_,low),(_,high) in itertools.product(primary[primary.frequency_Hz==.2].iterrows(),primary[primary.frequency_Hz==.5].iterrows()):
        fit,_=fit_gains([.2,.5],[low.gain_acceleration_ratio,high.gain_acceleration_ratio])
        fits.append(dict(output='head',start_s=5,end_s=35,selection=f'{low.trial} + {high.trial}',gain_0p5_over_0p2=high.gain_acceleration_ratio/low.gain_acceleration_ratio,**fit))
    pd.DataFrame(fits).to_csv(OUT/'head_first_order_fits.csv',index=False)
    pd.DataFrame(profiles).to_csv(OUT/'head_T_profile.csv',index=False)
    for f,block in primary.groupby('frequency_Hz'):
        axes[0].scatter(block.frequency_Hz,block.gain_acceleration_ratio,color='C0')
    dense=np.linspace(.12,.65,200)
    best=fits[0]
    if best['status']=='T_infinite_limit':
        pred=best['gamma_over_T_per_s']/(2*np.pi*dense)
        label='T -> infinity limit (no finite estimate)'
    else:
        pred=best['gamma_dimensionless']/np.sqrt(1+(2*np.pi*dense*best['T_s'])**2)
        label='First-order gain fit'
    axes[0].plot(dense,pred,label=label,color='C1')
    axes[0].set(xlabel='Frequency [Hz]',ylabel='Head / visual acceleration amplitude',title='Head output only: 5-35 s')
    axes[0].legend(fontsize=8)
    axes[1].set(xlabel='T [s], gamma optimized at each T',ylabel='SSE / minimum SSE',ylim=(.98,5),title='Dependence on analysis window')
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'head_gain_and_T_profile.png',dpi=180);plt.close(fig)
    report='''並進刺激0.2/0.5 Hzの同期不要ゲイン解析（2026-11-24）

入力: 両条件ともユーザー提示の変位片振幅0.0246 m。正弦波、固定床上の視覚並進刺激を仮定。
物理的な支持面移動ならモデルを変更する必要がある。
対象: gridtrans02/05、seguchi00001/00002の4試行。条件周波数はユーザー確認済みのファイル名に従う。
抽出済みAP頭部・関節モード加速度を使用。元CSVの記載先は現PCには存在せず、COM再構成用床反力は未入手。
モード和が頭部加速度と一致することを全6000点で確認したが、COMと同一ではない。

方法: 主区間5<=t<35秒（既存解析に合わせた区間であり、刺激ログで定常状態を確認したものではない）。
定数、直線トレンド、刺激周波数のcos/sinを同時最小二乗し、片振幅sqrt(cos係数^2+sin係数^2)を抽出。
mm/s²をm/s²に変換し、0.0246*(2*pi*f)^2で除算。位相は使用しない。
harmonic_partial_R2は直線トレンドのみの残差に対する刺激周波数成分追加の説明率。SNRやコヒーレンスではない。
別の5区間でも同じ解析を実施。区間比較は信頼区間ではない。

参考モデル: 頭部加速度/視覚加速度 = -gamma/(1+Ts)。T>=0、gamma>=0、gammaは無次元。
4試行の線形ゲイン二乗誤差を同じ重みで最小化。T=0およびT→∞の極限を明示的に比較。
COM・姿勢制御器を含む既存モデルのgamma（秒）とは定義も単位も異なり、比較・代用不可。
このモデルの0.5/0.2 Hzゲイン比は0.4〜1。2周波数平均の比がこの範囲外なら有限のT,gammaによる一致は不可能。
4通りの試行組合せは対応づけを仮定しない感度解析。2点の完全一致はモデルの妥当性の証拠にはならない。
繰り返し2試行のみで母集団の信頼区間は報告しない。
一定の時計オフセットには振幅が不変だが、時計ドリフト、非定常性、刺激振幅の実際の実装誤差には影響される。

主解析（頭部）の結果:
'''
    report+=primary[['trial','frequency_Hz','output_acceleration_peak_m_s2','gain_acceleration_ratio','harmonic_partial_R2']].to_string(index=False)
    report+='\n\n'+pd.DataFrame(fits).to_string(index=False)
    report+='\n\n結論: 頭部の参考解析のみ。従来のCOMモデルによるT,gammaの同定は床反力元CSVの入手待ち。\n'
    (OUT/'analysis_notes.txt').write_text(report,encoding='utf-8-sig')
    print(primary[['trial','frequency_Hz','gain_acceleration_ratio','harmonic_partial_R2']].to_string(index=False))
    print(pd.DataFrame(fits).to_string(index=False))
    print('Saved:',OUT)


if __name__=='__main__':
    main()
