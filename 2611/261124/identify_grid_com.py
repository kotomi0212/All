"""Unsynchronized COM gain identification; numpy/pandas only.

Run from any directory: python identify_grid_com.py --source PATH
The fixed-floor visual-stimulus assumption and physical calibration remain provisional.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import itertools
import numpy as np
import pandas as pd

G = 9.80665
A = .0246
WINDOWS = [(5,35),(5,25),(10,30),(15,35),(10,50),(0,60)]
OUT = Path(__file__).resolve().parent/'grid_com_identification'


def read_source(path):
    raw = path.read_bytes()
    rows = list(csv.reader(io.StringIO(raw.decode('cp932'))))
    indices = [i for i,r in enumerate(rows) if 'rawFx' in r]
    i = indices[-1]
    header = rows[i].copy(); header[0] = 'time_s'
    units = rows[indices[0]+1]
    assert len(set(header)) == len(header)
    d = pd.DataFrame([[float(v) for v in r] for r in rows[i+1:] if r], columns=header)
    assert np.isfinite(d).all().all() and np.allclose(np.diff(d.time_s), .01)
    meta = {r[0]:r[1] for r in rows[:indices[0]] if len(r)==2}
    assert float(meta['記録周波数(Hz)']) == 100
    assert units[header.index('Fx')] == 'N' and units[header.index('My')] in ('N・m','N\uff65m')
    assert meta['性別'] == '女性'
    return d, dict(path=str(path.resolve()),sha256=hashlib.sha256(raw).hexdigest(),
                  metadata=meta,units=dict(zip(header,units)),samples=len(d))


def harmonic(t,y,f):
    design = np.column_stack([np.ones(len(t)),t-t.mean(),np.cos(2*np.pi*f*t),np.sin(2*np.pi*f*t)])
    beta = np.linalg.lstsq(design,y,rcond=None)[0]
    residual = y-design@beta
    trend = design[:,:2]@np.linalg.lstsq(design[:,:2],y,rcond=None)[0]
    return beta[2]-1j*beta[3], float(1-np.sum(residual**2)/np.sum((y-trend)**2))


def plant(f,H,M):
    """Magnitude of existing VR velocity -> body COM displacement model, gamma=1,T=0."""
    mb=.978*M; jb=.0481*M*H**2; lb=.504*H
    inertia=jb/lb+mb*lb; kp=10.53*inertia; kd=3.97*inertia
    s=2j*np.pi*np.asarray(f); delay=np.exp(-.155*s)
    return abs(-kp*delay/((kp+kd*s)*delay-mb*G))


def fit(f,y,q):
    """Exact two-frequency solution plus T=0 and infinity; equal trial weights.

    Within-frequency repeated trials have identical plant magnitudes. Their means
    determine the optimum two-frequency prediction; their scatter remains in SSE.
    """
    f,y,q=map(np.asarray,(f,y,q)); w=2*np.pi*f
    levels=np.unique(f); assert len(levels)==2
    qq=[q[f==v] for v in levels]
    assert all(np.allclose(a,a[0]) for a in qq)
    def conditional(t):
        b=q/np.sqrt(1+(w*t)**2)
        gamma=float(b@y/(b@b))
        return float(np.sum((gamma*b-y)**2)),gamma
    r=(y[f==levels[1]].mean()/y[f==levels[0]].mean())/(qq[1][0]/qq[0][0])
    lo=levels[0]/levels[1]
    if r>=1:
        t=0.; sse,gamma=conditional(t); status='T_zero_limit'; rate=np.nan
    elif r<=lo:
        b=q/w; rate=float(b@y/(b@b)); sse=float(np.sum((rate*b-y)**2))
        t=gamma=np.inf; status='T_infinite_limit'
    else:
        t=float(np.sqrt((1-r*r)/(r*r*(2*np.pi*levels[1])**2-(2*np.pi*levels[0])**2)))
        sse,gamma=conditional(t); rate=np.nan; status='finite'
    return dict(T_s=t,gamma_s=gamma,gamma_over_T=rate,status=status,SSE=sse,
                normalized_gain_ratio=r,feasible_ratio_lower=lo,feasible_ratio_upper=1.),conditional


def verify():
    t=np.arange(3000)/100
    h,r=harmonic(t,.1+.002*t+.007*np.cos(2*np.pi*.2*t+.83),.2)
    assert abs(abs(h)-.007)<1e-12 and abs(r-1)<1e-12
    f=np.array([.2,.2,.5,.5]); q=plant(f,1.48,43.4)
    for tt in [.02,1.2,30]:
        y=.04*q/np.sqrt(1+(2*np.pi*f*tt)**2)
        result,_=fit(f,y,q)
        assert np.isclose(result['T_s'],tt) and np.isclose(result['gamma_s'],.04)
    assert fit(f,np.ones(4)*q,q)[0]['status']=='T_zero_limit'
    assert fit(f,q*np.where(f==.2,1,.1),q)[0]['status']=='T_infinite_limit'


def main(source):
    verify(); OUT.mkdir(exist_ok=True)
    amplitudes=[]; audits=[]; checks=[]
    paths=sorted(source.glob('*gridtrans0[25]_seguchi0000[12].csv'))
    assert len(paths)==4
    for path in paths:
        d,audit=read_source(path); audits.append(audit)
        meta=audit['metadata']; H=float(meta['身長(cm)'])/100; M=float(meta['体重(kg)']); mb=.978*M
        jb=.0481*M*H**2; lb=.504*H; c=jb/(mb*lb)+.039*H+lb
        f=.2 if 'trans02' in path.name else .5
        time=d.time_s.to_numpy(); rx=d.Fx.to_numpy(); ny=d.My.to_numpy()
        signals={'eq3_Fx_My':(ny-c*rx)/(mb*G),
                 'eq3_origin_minus12p3mm':(ny-(c+.0123)*rx)/(mb*G),
                 'eq3_origin_plus12p3mm':(ny-(c-.0123)*rx)/(mb*G),
                 'eq3_opposite_My':(-ny-c*rx)/(mb*G),
                 'vendor_xb':d.xb.to_numpy()*.001,
                 'force_equivalent_displacement':None}
        series={'time_s':time,'Fx_N':rx,'My_Nm':ny,'force_body_acceleration_m_s2':-rx/mb}
        series.update({k+'_m':v for k,v in signals.items() if v is not None})
        pd.DataFrame(series).to_csv(OUT/(path.stem+'_com.csv'),index=False)
        for start,end in WINDOWS:
            keep=(time>=start)&(time<end); t=time[keep]
            ah,ar=harmonic(t,-rx[keep]/mb,f)
            for name,values in signals.items():
                if values is None:
                    xh=ah/(-(2*np.pi*f)**2); r2=ar
                else:
                    xh,r2=harmonic(t,values[keep],f)
                amplitudes.append(dict(trial=path.stem,output=name,frequency_Hz=f,start_s=start,end_s=end,
                    samples=len(t),height_m=H,mass_kg=M,plant_gain=plant([f],H,M)[0],
                    displacement_peak_m=abs(xh),input_velocity_peak_m_s=A*2*np.pi*f,
                    measured_gain_s=abs(xh)/(A*2*np.pi*f),harmonic_partial_R2=r2,
                    force_acceleration_peak_m_s2=abs(ah),
                    reconstructed_over_force_acceleration_gain=abs(-(2*np.pi*f)**2*xh/ah),
                    reconstructed_over_force_phase_deg=np.angle(-(2*np.pi*f)**2*xh/ah,deg=True)))
        keep=(time>=5)&(time<35)
        for name,ours,vendor in [('acceleration',-rx/mb,d.axb.to_numpy()*.001),
                                 ('displacement',signals['eq3_Fx_My'],signals['vendor_xb']),
                                 ('cop',ny/(M*G),d.xp.to_numpy()*.001)]:
            x=ours[keep]-ours[keep].mean(); y=vendor[keep]-vendor[keep].mean()
            checks.append(dict(trial=path.stem,comparison=name,correlation=np.corrcoef(x,y)[0,1],
                               demeaned_rmse=float(np.sqrt(np.mean((x-y)**2))),
                               vendor_on_candidate_slope=float(x@y/(x@x))))
        audit['vertical_load_mass_kg']=float(d.Fz.mean()/G)
        audit['raw_vs_processed_horizontal_mean_N']=[float(d.rawFx.mean()),float(d.Fx.mean())]
    data=pd.DataFrame(amplitudes); data.to_csv(OUT/'com_harmonic_amplitudes.csv',index=False)
    fits=[]; profiles=[]
    for (name,start,end),block in data.groupby(['output','start_s','end_s'],sort=False):
        selections=[('all_four_trials',block)]
        if (start,end)==(5,35):
            for (_,low),(_,high) in itertools.product(block[block.frequency_Hz==.2].iterrows(),block[block.frequency_Hz==.5].iterrows()):
                selections.append((low.trial+' + '+high.trial,pd.DataFrame([low,high])))
        for selection,b in selections:
            result,conditional=fit(b.frequency_Hz,b.measured_gain_s,b.plant_gain)
            fits.append(dict(output=name,start_s=start,end_s=end,selection=selection,**result))
            if selection=='all_four_trials':
                for tv in np.logspace(-4,5,400):
                    ss,gg=conditional(tv)
                    profiles.append(dict(output=name,start_s=start,end_s=end,T_s=tv,gamma_s=gg,SSE=ss,
                                         SSE_over_minimum=ss/result['SSE'] if result['SSE']>0 else np.nan))
    fits=pd.DataFrame(fits); fits.to_csv(OUT/'com_parameter_fits.csv',index=False)
    pd.DataFrame(profiles).to_csv(OUT/'com_T_profile.csv',index=False)
    pd.DataFrame(checks).to_csv(OUT/'vendor_consistency.csv',index=False)
    (OUT/'source_audit.json').write_text(json.dumps(audits,ensure_ascii=False,indent=2),encoding='utf-8')
    primary=data[(data.output=='eq3_Fx_My')&(data.start_s==5)&(data.end_s==35)]
    selected=fits[(fits.output=='eq3_Fx_My')&(fits.selection=='all_four_trials')]
    note='''元CSVによるCOM T・γゲイン同定（2026-11-24）

結論: 既存のVR速度→身体部COM変位モデルで、安定した有限のT・γは確定できない。
主区間5<=t<35秒の4試行同時推定はT→∞、γ→∞の極限。
gamma_over_Tはその極限での比であり、個別のT・γの推定値ではない。
探索上限を同定値として報告していない。区間別に有限値も現れるが再現しない。

入力・元データ:
Desktop/1022csvのgridtrans02/05、各00001/00002の4ファイル。
各6000点、100 Hz、0..59.99秒、cp932。SHA256と全ヘッダーはsource_audit.json。
身長148cm、女性、体重は試行順に43.4/43.5/43.5/43.7kg（各CSV値を使用）。
Fz平均/gは約44.45/44.50/44.53/44.66kgで、メタデータ体重とは一致しない。
ファイル名による0.2/0.5Hz、視覚変位片振幅0.0246mを引き継ぐ。
内部試験名は4件ともgridtrans02。入力の位相は未同期のため使用しない。
固定床上の視覚並進を暫定仮定。実際に床が動く場合はこのモデルを採用できない。
刺激開始終了ログはCSVにない。5..35秒は既存解析と比較する主区間であり、
刺激定常区間を独立確認したものではない。他の5区間は感度解析。

COM再構成と同定:
処理済みFx/My列の単位N/N･mを確認。AP候補Rx=Fx,Ny=My。
My/(M*g)とCSVのAP COP=xpの相関は0.9988以上で、この軸候補と整合。
ただしこれは物理的な符号校正の証明ではない。
mb=.978M,Jb=.0481MH²,lb=.504H,Lf=.039H（引継ぎの2019女性係数）。
xb=(My-c*Fx)/(mb*g), c=Jb/(mb*lb)+Lf+lb。gを含める。
CSVのxb列は比較用。別の2025身体係数へ暗黙に切り替えていない。
モーメント原点12.3mmの補正有無が確定しないため0と±12.3mmを保存。
さらにMy反転も感度解析として保存し、適合が良いことを根拠に符号を選ばない。
Fx/Myはエクスポート済み列。raw接頭辞列も校正前の生信号とは断定しない。
CSVヘッダーにHP0.1Hz/LP4Hzの記載あり。列別の適用処理・次数は不明。
追加フィルタと未知のフィルタ応答の逆補正は適用していない。

定数+直線トレンド+cos+sinの同時最小二乗で刺激周波数の片振幅を抽出。
出力COM変位振幅/(0.0246*2πf)が入力速度→COM変位ゲイン[s]。
Gx=-kp*γ*exp(-td*s)/[(1+T*s)*{(kp+kd*s)*exp(-td*s)-mb*g}]。
kp=10.53*(Jb/lb+mb*lb),kd=3.97*(Jb/lb+mb*lb),td=.155秒。
従来どおり制御モデルの慣性項を省略。Tとγの単位はともに秒。
4試行の線形ゲインSSEを等重みで最小化。二周波数の解析解とT=0/∞を比較。
同じ身長・係数では単位γ,T=0のプラント振幅qは質量によらず一致する。
正規化ゲイン比=(平均gain_0.5/平均gain_0.2)/(q_0.5/q_0.2)。
有限T>0で可能な比は0.4より大きく1未満。主区間の比0.338491は範囲外。
harmonic_partial_R2はトレンド残差に対する刺激周波数追加の説明率。
信頼区間・コヒーレンス・時系列モデル適合率ではない。

独立の力学的一致チェック:
足部静止条件で身体部COM加速度=-Fx/mbを別に算出。
その周波数係数/[-(2πf)²]はforce_equivalent_displacementとして比較。
これは加速度の周波数領域換算であり、直接測定変位ではない。
主区間のEq3変位から換算した加速度/床反力加速度の振幅比は0.038..0.212、
位相差も約-107,+121,+123,+14度で一致しない。
同一CSV内の相対位相なので刺激との時計同期は不要。
この不一致の原因は、列別処理、軸符号、単一剛体近似、非定常性などを分離できない。
DC原点や身長体重を調整して整合性を保証したものではない。
My反転では主区間T約0.428秒、γ約0.0145秒になるが採用根拠にはしない。

次に必要な証拠:
固定床の視覚刺激であること、刺激の開始終了（記録時間との大まかな対応で可）、
当日の被験者方向、既知荷重によるFx/My符号・原点、列別フィルタ処理。
これらが得られれば区間・軸・処理を確定し再実行する。

検証: 未知位相/オフセット/ドリフトを含む合成波の振幅回復、
既知T=.02/1.2/30秒・γ=.04秒の回復、T=0/∞の境界判定を実行時assertで確認。
元CSV、既存頭部解析、旧結果は変更していない。

主区間の振幅・ゲイン:
'''
    note+=primary[['trial','displacement_peak_m','measured_gain_s','harmonic_partial_R2']].to_string(index=False)
    note+='\n\n区間別4試行同時推定:\n'+selected.to_string(index=False)
    note+='\n\n主区間の出力定義・符号・原点感度:\n'+fits[(fits.start_s==5)&(fits.end_s==35)&(fits.selection=='all_four_trials')].to_string(index=False)
    (OUT/'COM_ANALYSIS.txt').write_text(note+'\n',encoding='utf-8-sig')
    print(primary[['trial','measured_gain_s','harmonic_partial_R2']].to_string(index=False))
    print(selected.to_string(index=False))
    print(pd.DataFrame(checks).to_string(index=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--source',type=Path,default=Path('C:/Users/seguchi/Desktop/1022csv'))
    main(parser.parse_args().source)
