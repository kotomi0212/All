"""Reanalyze explicitly selected ten-tone records, keeping trials/conditions separate."""
import argparse,csv,json,subprocess,sys,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_periodic_response import write

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('trials',nargs='+',type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    subjects={};manifest=[];allfrf=[];allcycles=[];inputs=[];summaries=[];before={}
    for index,folder in enumerate(args.trials):
        folder=folder.resolve();meta=json.loads((folder/(folder.name+'_metadata.json')).read_text(encoding='utf-8-sig'))
        subject=meta['participant']['subject']
        if subject not in subjects:subjects[subject]='Participant_'+chr(65+len(subjects))
        alias=subjects[subject];label=alias+'_T'+str(meta['participant']['trial_number'])
        condition=meta['stimulus']['condition'];tag=f'{index+1:02d}_{label}'
        for src in folder.iterdir():
            if src.is_file():before[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest()
        dest=out/tag
        subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).with_name('analyze_periodic_response.py')),'--trials',str(folder),'--label',label,'--output',str(dest)],check=True)
        summary=json.loads((dest/'summary.json').read_text());summaries.append(summary['subjects'][0])
        manifest.append(dict(participant=alias,trial_number=meta['participant']['trial_number'],trial_id=meta['trial_id'],condition=condition,source=str(folder),output=tag,termination=meta['termination'],analysis_start_s=10,analysis_end_s=90,period_s=20,cycles=4,imu_status=meta.get('external_imu_quality',{}).get('status'),mass_kg=meta.get('weight_estimate',{}).get('body_mass_kg')))
        for name,target in [('frequency_response.csv',allfrf),('cycle_coefficients.csv',allcycles),('input_log_checks.csv',inputs)]:
            for row in read(dest/name):target.append(dict(**row,condition=condition,recording_tag=tag,subject_alias=alias))
    write(out/'trial_manifest.csv',manifest);write(out/'frequency_response_all.csv',allfrf);write(out/'cycle_coefficients_all.csv',allcycles);write(out/'input_checks_all.csv',inputs)
    fig,axs=plt.subplots(len(manifest),2,figsize=(12,3*len(manifest)),squeeze=False,layout='constrained')
    lines=['# 100秒記録の再解析','',f'対象{len(manifest)}試行。全試行で10–90秒、20秒×4周期。既存の時刻対応を使用し、相関最大化シフトなし。条件は混合せず各試行を個別解析。','', '|被験者|試行|刺激条件|頭部IMU|','|---|---|---|---|']
    for row in manifest:lines.append(f"|{row['participant']}|{row['trial_number']}|{row['condition']}|{row['imu_status']}|")
    lines+=['','速度振幅一定と加速度振幅一定では入力の周波数構成の強さが異なる。生振幅やS/Nを条件の違いを無視して人物差として比較しない。','', '|被験者・試行|条件|出力|S/N中央値 dB|コヒーレンス中央値|再現性が相対的に高い点|','|---|---|---|---:|---:|---|']
    for index,(m,s) in enumerate(zip(manifest,summaries)):
        tag=m['output'];cond='equal velocity' if 'equal_velocity' in m['condition'] else 'equal acceleration'
        for name,color,label in [('COM_acceleration_Fy_over_M','#d62728','COM acceleration'),('Head_Z_raw','#0072B2','Head IMU Z (raw)')]:
            rr=[r for r in allfrf if r['recording_tag']==tag and r['output']==name and r.get('gain')]
            if not rr:continue
            f=np.array([float(r['frequency_hz']) for r in rr]);coh=np.array([float(r['coherence']) for r in rr]);snr=np.array([float(r['snr_db']) for r in rr])
            axs[index,0].semilogx(f,coh,'o-',color=color,label=label);axs[index,1].semilogx(f,snr,'o-',color=color)
            best=sorted(rr,key=lambda r:float(r['coherence']),reverse=True)[0]
            lines.append(f"|{s['participant']}|{m['condition']}|{label}|{np.median(snr):.2f}|{np.median(coh):.3f}|{best['frequency_hz']} Hz：C={float(best['coherence']):.2f}, S/N={float(best['snr_db']):.1f} dB|")
        for j,ax in enumerate(axs[index]):
            ax.set_title(s['participant']+' | '+cond);ax.grid(alpha=.2,which='both');ax.set_xlabel('Frequency (Hz)');ax.set_ylabel('Coherence (4 cycles)' if j==0 else 'S/N (dB)')
        axs[index,0].set_ylim(0,1);axs[index,0].legend(fontsize=7)
    fig.suptitle('100 s records | each trial/condition analyzed separately | gray cycle values in individual FRF figures')
    fig.savefig(out/'all_trials_reproducibility.png',dpi=150);plt.close(fig)
    repeat=[]
    for subject in subjects.values():
        for condition in {m['condition'] for m in manifest}:
            group=[m for m in manifest if m['participant']==subject and m['condition']==condition]
            if len(group)<2:continue
            for i,one in enumerate(group[:-1]):
                for two in group[i+1:]:
                    a=[r for r in allfrf if r['recording_tag']==one['output'] and r['output']=='COM_acceleration_Fy_over_M']
                    b=[r for r in allfrf if r['recording_tag']==two['output'] and r['output']=='COM_acceleration_Fy_over_M']
                    for x,y in zip(a,b):
                        h1=complex(float(x['H_real']),float(x['H_imag']));h2=complex(float(y['H_real']),float(y['H_imag']))
                        repeat.append(dict(participant=subject,condition=condition,trial1=one['trial_number'],trial2=two['trial_number'],frequency_hz=x['frequency_hz'],gain_ratio_trial2_over_trial1=abs(h2/h1),phase_difference_deg=np.angle(h2/h1,deg=True),coherence_trial1=x['coherence'],coherence_trial2=y['coherence']))
    write(out/'between_trial_comparison.csv',repeat)
    lines+=['','## 解釈と限界','', 'S/Nは各周期Z=Y/Uの複素平均パワー/周期間の複素残差分散を10log10で換算したもの。4周期のコヒーレンスと円周位相ばらつきもCSVに記録。複素平均の有限反復バイアス補正版を併記。単一周期からコヒーレンスを出していない。','全ゼロIMUは応答ゼロとは判断せず、記録なしとして評価不可。利用可能な頭部IMUも重力補正の妥当性が未確認のためZ軸生カウントの補助出力。COM加速度=Fy/M（符号未校正、Mは静止Fz推定）。COM位置は必要な身体パラメータ未確定のため仮値で再構成せず評価不可。','複数試行のある同一被験者・同一条件だけで試行間ゲイン比と位相差を別CSVに出した。独立試行2回では再現性の確定や信頼区間は難しく、周期間のK=4と試行間のn=2を混同しない。','kp,kd,tdの個別推定と位相符号・物理遅延の確認が不足するため、原著のTv/γvのフィットは行わない。既存モデルの一律定数を個人の推定結果として扱わない。','応答の高い点が点在していても連続した同定可能帯域とは断定しない。試行別のゲイン・位相・周期散布・S/N・コヒーレンス図で確認する。','', '入力ログの10周波数・振幅・微分関係、欠損・時刻情報は各試行のinput_log_checks.csv、synchronization_audit.csv、signal_schema.csvに保存。元ファイルは変更せずSHA256を照合。','', '## 再実行','', '```powershell', '.\\.venv\\Scripts\\python.exe -X utf8 research_sync/analysis_tools/reanalyze_recordings.py --output research_sync/results/NEW_OUTPUT TRIAL_FOLDER_1 TRIAL_FOLDER_2 ...', '```']
    (out/'RESULTS.md').write_text('\n'.join(lines),encoding='utf-8')
    for p,digest in before.items():
        if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=digest:raise ValueError('Source changed')
    status=dict(trials=len(manifest),cycles=4*len(manifest),frequency_rows=len(allfrf),cycle_rows=len(allcycles),source_files_unchanged=len(before),between_trial_rows=len(repeat))
    (out/'completion_check.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))

if __name__=='__main__':main()
