"""Historical rest reference and a Japanese report for cycle FRF outputs."""
import argparse,json,csv,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_periodic_response import table,harmonic,write

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    root=args.root.resolve();out=args.output.resolve()
    summary=json.loads((out/'summary.json').read_text(encoding='utf-8'))
    with (out/'trial_inventory.csv').open(encoding='utf-8-sig') as f:inventory=list(csv.DictReader(f))
    with (out/'frequency_response.csv').open(encoding='utf-8-sig') as f:frf=list(csv.DictReader(f))
    with (out/'input_log_checks.csv').open(encoding='utf-8-sig') as f:checks=list(csv.DictReader(f))
    with (out/'source_hashes.csv').open(encoding='utf-8-sig') as f:hashes=list(csv.DictReader(f))
    for r in hashes:
        if hashlib.sha256((root/r['file']).read_bytes()).hexdigest()!=r['sha256']:raise ValueError('Original source changed')
    F=np.array(summary['subjects'][0]['frequency_hz']);baseline=[];baseinv=[]
    fig,axs=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    rest=[r for r in inventory if 'rest' in r['condition'] and 'duplicate' not in r['reason']]
    for index,item in enumerate(rest):
        mp=root/item['metadata'];m=json.loads(mp.read_text(encoding='utf-8-sig'));folder=mp.parent;stem=mp.name.removesuffix('_metadata.json')
        note=dict(reference='Rest '+str(index+1),trial_id=m['trial_id'],status='unavailable',reason='',matched_main_80s='unavailable: rest duration 60 s',comparison='20 s matches main cycle only; day, orientation and movements not controlled')
        try:
            if m.get('termination')!='complete':raise ValueError('Incomplete rest trial')
            q,_=table(folder/(stem+'_quest.csv'));b,_=table(folder/(stem+'_combalance.csv'))
            t=np.interp(b['pc_estimated_sample_time_s'],q['pc_mapped_time_s'],q['stimulus_elapsed_s'],left=np.nan,right=np.nan)
            mass=(m.get('weight_estimate') or {}).get('body_mass_kg')
            for start in [10.,30.]:
                k=(t>=start)&(t<start+20)
                if k.sum()<1900:continue
                signals={'Head_Z_raw':(b['external_imu_count_2'],'counts',1)}
                if mass and mass>0:signals['COM_acceleration_Fy_over_M']=(b['sensor_Fy_N']/mass,'m/s^2',0)
                for name,(y,unit,panel) in signals.items():
                    a=abs(harmonic(t[k],y[k],F))
                    axs[panel].loglog(F,a,'o-',alpha=.5,label=f'Rest {index+1}, {start:g}–{start+20:g} s')
                    for hz,v in zip(F,a):baseline.append(dict(reference=note['reference'],trial_id=m['trial_id'],output=name,start_s=start,end_s=start+20,frequency_hz=hz,amplitude=v,unit=unit,status='unmatched historical reference only'))
            note.update(status='reference only',reason='No 80 s baseline; cannot assign observed fluctuations to spontaneous sway vs sensor noise or IMU checks')
        except Exception as e:note['reason']=str(e)
        baseinv.append(note)
    for ax,title,unit in zip(axs,['COM acceleration reference','Head Z raw reference'],['m/s²','counts']):
        ax.set(title=title,xlabel='Frequency (Hz)',ylabel='Amplitude ('+unit+')');ax.grid(alpha=.25);ax.legend(fontsize=6)
    fig.suptitle('Unstimulated historical records: 20 s windows; unmatched orientation/day; not S/N denominator')
    fig.savefig(out/'rest_reference.png',dpi=150);plt.close(fig)
    write(out/'baseline_inventory.csv',baseinv);write(out/'baseline_amplitudes.csv',baseline)
    comparisons=[dict(participant=r['participant'],frequency_hz=r['frequency_hz'],output='COM_position_Sonobe',status='unavailable',reason=r['reason']) for r in frf if r['output']=='COM_position_Sonobe']
    write(out/'output_comparison.csv',comparisons)
    fig,axs=plt.subplots(3,2,figsize=(12,10),sharex=True,layout='constrained')
    for index,subject in enumerate(sorted(summary['subjects'],key=lambda x:x['participant'])):
        for name,color,label in [('COM_acceleration_Fy_over_M','#d62728','COM acceleration'),('Head_Z_raw','#0072B2','Head IMU Z (raw)')]:
            rr=[r for r in frf if r['participant']==subject['participant'] and r['output']==name]
            hz=[float(r['frequency_hz']) for r in rr]
            axs[index,0].semilogx(hz,[float(r['coherence']) for r in rr],'o-',color=color,label=label)
            axs[index,1].semilogx(hz,[float(r['snr_db']) for r in rr],'o-',color=color,label=label)
        axs[index,0].set(title=subject['participant'],ylabel='Coherence (4 cycles)',ylim=(0,1))
        axs[index,1].set(title=subject['participant'],ylabel='S/N (dB; cycle residual)')
        for ax in axs[index]:ax.grid(alpha=.2,which='both');ax.set_xticks(F,labels=[str(x) for x in F],rotation=60,fontsize=8)
    axs[0,0].legend(fontsize=8);axs[-1,0].set_xlabel('Frequency (Hz)');axs[-1,1].set_xlabel('Frequency (Hz)')
    fig.suptitle('Within-trial reproducibility | 1 trial per participant, 4 complete cycles | no pass/fail threshold')
    fig.savefig(out/'reproducibility_overview.png',dpi=150);plt.close(fig)
    lines=['# 速度振幅一定10周波数：周期間の応答再現性','',
      '## 計算できた結果','',
      f"主条件の完了試行は{summary['analyzed_trials']}名・各1試行。各10–90秒、20秒周期×4回、合計12周期。試行間再現性は評価不可。",
      '入力周波数は設定とログのFFTピークで一致：'+', '.join(map(str,F))+' Hz。FFTの丸め値で回帰せず、設定の正確な周波数と実タイムスタンプを使用。',
      f"入力速度振幅の設定値からの最大相対差は{100*max(float(r['relative_error_vs_config']) for r in checks):.3f}%。",
      '', '|被験者|出力|全10点のS/N中央値 dB|コヒーレンス中央値|探索的BH q<0.05の点 Hz|','|---|---|---:|---:|---|']
    for s in sorted(summary['subjects'],key=lambda x:x['participant']):
        for o in s['outputs']:
            lines.append(f"|{s['participant']}|{o['output']}|{o['median_snr_db']:.2f}|{o['median_coherence']:.3f}|{o['q_below_0p05_hz']}|")
    lines+=['','S/N中央値は周波数全体を要約した記述値で、被験者間の平均ではない。コヒーレンスのp/qは「4周期が独立な複素ガウス雑音」という近似の探索的値。連続周期の依存、周波数間の依存、小標本を考慮すると確証ではない。一律の合格閾値にはしない。','',
      '## 被験者ごとの解釈','']
    for s in sorted(summary['subjects'],key=lambda x:x['participant']):
        subset=[r for r in frf if r['participant']==s['participant'] and r['output']=='COM_acceleration_Fy_over_M']
        strongest=sorted(subset,key=lambda r:float(r['coherence']),reverse=True)[:3]
        lines.append(s['participant']+'：床反力COM加速度で相対的に再現性が高い点は '+', '.join(f"{r['frequency_hz']} Hz (C={float(r['coherence']):.2f}, S/N={float(r['snr_db']):.1f} dB, 位相SD={float(r['phase_circular_sd_deg']):.1f}°)" for r in strongest)+'。点在する候補であり、連続した帯域で再現よく検出できたとは結論しない。')
    lines+=['','床反力COM加速度を次の主出力候補とする。F=Maに対応し、位置のモデル誤差や頭部IMUの重力補正に依存しないため。ただし本データで十分なS/Nが確認されたという意味ではない。COM位置は個別身体パラメータ未確定のため比較不可。出力比較CSVに理由を記録。位置を後日推定しても床反力を共有するため独立な精度検証にはならない。','',
      '## 座標・単位・同期の前提','',
      '`quest.js` は world.position.z=base_z-s.x として背景を動かし、s.vをログへ保存。前方=-world Zの規約なら背景前進速度はログの+ap_velocity_m_s。仮想カメラの移動ではないため追加の符号反転はしない。頭部の実位置ログではない。',
      '床反力はsensor_Fy_N/M、Mは試行前静止Fzの絶対値中央値/9.80665で推定（装着物を含み、独立校正なし）。+Fyが身体前方への合力とみなせることは実機で未検証。振幅は符号反転に不変だが、位相には180°の不確定性。身体の前向きとプレート座標の実機照合が必要。',
      'Quest pc_mapped_time_s と床反力 pc_estimated_sample_time_s をそのまま用い、補間で刺激経過時間に対応づけた。出力の相関最大化シフトなし。クロック残差と欠損は同期監査CSV参照。PC受信に基づく床反力時刻であり、センサー取得・表示発光遅延の校正を意味しない。',
      '頭部IMUはexternal_imu_count_2（ユーザー指定のセンサーZ=AP）を生カウントの補助出力とした。全7chは列一覧に収録。3軸加速度・角速度のデコードは既存コードとメーカーAPI設定を参照できるが、姿勢補正の信頼性が未確認のため、実験室座標頭部加速度や頭部変位とは扱わない。',
      'COM位置の確定値は評価不可：性別・個別身体パラメータ・モーメント原点・符号が未確定。以前の図では既定の男性/女性係数を仮定していたが、今回は要求に従いどちらも本人の値として採用せず、変位FRF・S/Nを算出しない。式は x=(-Mx-c*Fy)/(mb*g), c=Jb/(mb*lb)+Lf+lb, mb=.978M。男性 Jb=.0473MH², lb=.513H、女性 .0481MH², .504H、Lf=.039H。各感度結果は推定法依存。',
      '', '## 複素応答・S/N・コヒーレンスの定義','',
      '各20秒区間で、切片・直線ドリフト・10周波数すべてのcos/sinを同時最小二乗回帰。a cos + b sin の複素係数は a-jb。区間ごとの絶対刺激経過時間を使うので4周期の入力位相がそろう。窓関数・出力フィルターなし。最長周期20秒では1周期しかなく、直線ドリフト項との区別に限界がある。',
      '各周期Z_k=Y_k/U_k。代表FRFは H=sum(conj(U_k)Y_k)/sum(|U_k|²)。位相角を平均せず複素係数を集約。灰色は各周期、青は集約。',
      '複素平均Zbar、周期間分散 V=sum|Z_k-Zbar|²/(K-1)。S/N=10log10(|Zbar|²/V)。これは各周期における同期応答/残差の比で、平均後のS/Nではない。有限Kによる同期パワーの上振れを補正した値も併記：10log10((|Zbar|²-V/K)/V)。分子<=0は空欄（雑音から分離不可）。',
      'コヒーレンス=|sum(conj(U_k)Y_k)|²/[sum|U_k|² sum|Y_k|²]。K=4を明示。K=1は評価不可。位相ばらつきは単位複素数平均の長さRからsqrt(-2 ln R)を度へ換算。周期間の自然動揺・非定常な応答・計測ノイズを含み、それら個々の成分は分離しない。',
      '', '## 無刺激の参考','',
      f'静止条件{len(baseinv)}記録を確認。主解析と同じ80秒の無刺激記録はなく、20秒（主解析の1周期と同じ長さ）で参考値を算出した。日付・計測方向・IMU動作確認の混入を統制していないため、今回の3名の静止基準やS/N分母には使用しない。記録ごとの理由はbaseline_inventory.csv。',
      '', '## Sonobe & Miura (2025)への接続','',
      '式(8)の視覚経路は Tv*d(x_target)/dt+x_target=-gamma_v*v_background として保持する。固定床の支持面加速度は0。身体慣性を含む閉ループの一例は Gx=-kp*gamma_v*exp(-td*s)/[(Tv*s+1){A*s²+(kp+kd*s)exp(-td*s)-mb*g}], A=Jb/lb+mb*lb。加速度出力はGa=s²Gx。力の符号と入力規約を合わせる必要がある。既存の慣性項省略版や加速度一次遅れの直接フィットとは区別する。',
      '同一被験者について独立に同定したkp,kd,tdは確認できなかった。既存スクリプトの10.53,3.97,0.155 sは共通固定値で、個人の校正結果ではない。従ってTv/γvのフィット・Tv走査は今回は行わない。床反力の位相規約と物理遅延の確認も不足。一次遅れを実測COMへ直接当てはめた過去の結果を原著のTv/γvとして引用しない。',
      '原著情報：https://doi.org/10.1109/TNSRE.2025.3632867 。著者研究室の説明：https://ddlab-kut.net/archives/727 （0.05–0.70 Hzを対象、頭部加速度制御は対象外）。',
      '', '## 次に必要なデータ','',
      '0.05–0.70 Hzを優先して同条件の独立な試行を追加し、同じ日・同じ向きの無刺激試行を確保する。今回の4周期の結果では連続した同定帯域は確定できない。Participant Cの0.90 Hzなどの点状候補は低周波COMモデルの帯域外であり、それだけでTv同定へ進まない。',
      '', '## 実行と確認','',
      '```powershell',
      '.\\.venv\\Scripts\\python.exe -X utf8 research_sync/analysis_tools/analyze_periodic_response.py --output research_sync/results/NEW_OUTPUT',
      '.\\.venv\\Scripts\\python.exe -X utf8 research_sync/analysis_tools/periodic_reference_report.py --output research_sync/results/NEW_OUTPUT',
      '.\\.venv\\Scripts\\python.exe -m unittest discover -s research_sync/analysis_tools -p test_periodic_response.py',
      '```','',
      f'全対象{summary["analyzed_trials"]}完了試行で実行済み、周波数一覧 {len(frf)}行（計算60行・COM位置評価不可30行）、周期係数240行。入力・FRF・静止参考図を保存。元データのSHA256を再照合して変更なしを確認。中断試行、旧条件、取得コピーの除外理由をtrial_inventory.csvへ保存。数値不足で算出不可の値は空欄/NaNとし、0で埋めない。']
    (out/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({'report':str(out/'RESULTS.md'),'rest_records':len(baseinv),'baseline_rows':len(baseline),'source_hashes_verified':len(hashes)},indent=2))

if __name__=='__main__':main()
