"""Export unchanged device modal acceleration samples, with provenance checks."""
from pathlib import Path
import csv
import io
import html
import json
import zipfile
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from head_imu_analysis import BASE, NAMES

OUT = Path(__file__).resolve().parent/'outputs'/'joint_acceleration_20261101'
COLS = ['xix1','xix2','xiy1','xiy2']
REFERENCE = 'https://doi.org/10.3389/fneur.2025.1500308'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams['font.family'] = ['Yu Gothic','Meiryo','DejaVu Sans']
    checks, manifest, cards = [], [], []
    for name in sorted(NAMES):
        path = BASE/f'261022_{name}.csv'
        raw = path.read_bytes()
        try:
            text = raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            text = raw.decode('cp932')
        rows = list(csv.reader(io.StringIO(text)))
        idx = [i for i,r in enumerate(rows) if set(COLS).issubset(r)][-1]
        header = rows[idx]
        selected = [0]+[header.index(c) for c in COLS]
        samples = [r for r in rows[idx+1:] if r and any(v.strip() for v in r)]
        assert all(len(r) == len(header) for r in samples)
        # Copy original decimal strings: no rounding, centering or filtering.
        copied = [[r[j] for j in selected] for r in samples]
        data = np.asarray(copied, dtype=float)
        t = data[:,0]
        assert np.isfinite(data).all() and np.all(np.diff(t)>0)
        assert np.allclose(np.diff(t),.01,atol=1e-8)
        assert len(t)==6000
        unit_rows = [r for r in rows[:idx] if len(r)==len(header) and all(r[header.index(c)]=='mm/s/s' for c in COLS)]
        assert unit_rows
        dest=OUT/f'{path.stem}_joint_acceleration.csv'
        with dest.open('w',newline='',encoding='utf-8-sig') as f:
            w=csv.writer(f)
            w.writerow(['time_s','AP_ankle_xix1_mm_s2','AP_hip_xix2_mm_s2',
                        'ML_ankle_xiy1_mm_s2','ML_hip_xiy2_mm_s2'])
            w.writerows(copied)
        with dest.open(encoding='utf-8-sig',newline='') as f:
            exported=list(csv.reader(f))[1:]
        assert exported == copied
        for axis,a,h,hd in [('AP','xix1','xix2','axhd'),('ML','xiy1','xiy2','ayhd')]:
            x=np.array([[float(r[header.index(c)]) for c in [a,h,hd]] for r in samples])
            residual=x[:,0]+x[:,1]-x[:,2]
            checks.append([path.name,axis,len(t),float(np.sqrt(np.mean(residual**2))),
                           float(np.max(abs(residual))),float(np.std(x[:,2])),
                           float(np.std(residual)), 'mm/s^2'])
        meta={r[0]:r[1] for r in rows[:20] if len(r)==2}
        manifest.append(dict(source=str(path),sha256=hashlib.sha256(raw).hexdigest(),
                             internal_recording_name=rows[0][1],metadata=meta,
                             output=dest.name,samples=len(t),start_s=float(t[0]),last_s=float(t[-1]),
                             sampling_hz=100,additional_filter=False,mean_subtraction=False,
                             resampling=False,method='verbatim device-exported modal acceleration columns',
                             source_column_order=COLS,decimal_strings_verified=True))
        fig,axes=plt.subplots(2,2,figsize=(14,7),layout='constrained',sharex=True)
        fig.suptitle(f'{path.stem}\n足関節・股関節モードによる頭部加速度 | 装置出力値そのまま',fontsize=14)
        for j,(col,label,color) in enumerate(zip(COLS,
                ['AP 前後：足関節由来','AP 前後：股関節由来','ML 左右：足関節由来','ML 左右：股関節由来'],
                ['#167d9a','#d37828','#167d9a','#d37828'])):
            ax=axes.flat[j]
            ax.plot(t,data[:,j+1],lw=.7,color=color)
            ax.set(title=f'{label} ({col})',xlabel='時刻 (s)',ylabel='加速度 (mm/s²)',xlim=(0,60))
            ax.grid(alpha=.2)
        for pair in axes:
            lo=min(ax.get_ylim()[0] for ax in pair); hi=max(ax.get_ylim()[1] for ax in pair)
            for ax in pair: ax.set_ylim(lo,hi)
        fig.savefig(OUT/f'{path.stem}.png',dpi=140)
        plt.close(fig)
        cards.append(f'<h2>{html.escape(path.stem)}</h2><a href="{dest.name}">波形CSV</a><img src="{path.stem}.png">')
    with (OUT/'reconstruction_check.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f); w.writerow(['file','axis','samples','sum_minus_head_RMSE',
            'sum_minus_head_max_abs','head_demeaned_RMS','sum_minus_head_demeaned_RMS','unit']);w.writerows(checks)
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    notes='''足関節・股関節モードによる頭部加速度の時刻波形
参照: Sonobe & Mitsutake (2025), 2.2節、式(7)(8)、Table 1。
https://doi.org/10.3389/fneur.2025.1500308

前後(AP): xix1=足関節モード、xix2=股関節モード。
左右(ML): xiy1=足関節モード、xiy2=股関節（論文のモデルでは腰部）モード。
各値は各モードが生む頭部加速度（mm/s²）。関節点そのものの並進加速度やトルクではありません。
全13件・100 Hz・6000点、時刻0–59.99秒。全区間を出力。

今回の操作は装置が推定済みのモード列の抽出・可視化です。生IMU・床反力からモデルを再推定した結果ではありません。
CSV内の数値文字列をそのまま転記し、出力後に全78000行の文字列一致を確認。
追加フィルタなし、平均除去なし、ドリフト補正なし、再サンプリングなし、割合への正規化なし。
元CSVのメタデータにはHP 0.1 Hz / LP 4 Hzの設定があります。装置側の既存処理を解除した完全な生波形ではありません。
論文2.1節にも推定時のKalman filterの使用が記載されています。今回それを追加・再実装していません。

論文式(7): [AP足関節, AP股関節]^T = C_AP [下部COM加速度, 上部COM加速度]^T。
男性 C_AP=[[0.00917,1.37],[-2.06,0.903]]、女性 C_AP=[[0.0343,1.36],[-2.04,0.910]]。
論文式(8): [ML足関節, ML股関節]^T = C_ML [脚部COM加速度, 上部COM加速度]^T。
男性 C_ML=[[0.512,1.06],[-3.71,1.86]]、女性 C_ML=[[0.587,1.01],[-3.51,1.78]]。
上記は論文掲載の丸め係数を記録したもので、今回の抽出に掛け直していません。
CSVにはこれらの入力となる上下部COM加速度が明示的に保存されていないため、式(7)(8)による独立な再計算は行っていません。
列名と既存プロジェクトの対応付けに基づく装置モード波形であり、装置内部計算の完全検証ではありません。

整合性: 論文では両モードの和が頭部加速度に一致します。各記録の差をreconstruction_check.csvに記録。
APはほぼ一致しますがMLには差があります。差をモード側に分配したり補正したりしていません。
MLを論文どおりに独立再計算・検証するには装置の上下部加速度・身体パラメータ・座標変換と推定処理の仕様確認が必要です。
gridtrans05の2件は内部測定名がgridtrans02です。ファイル名を結果名とし、内部名もmanifestに保存。
原本は変更していません。manifest.jsonに原本SHA256、メタデータ、処理条件を保存。
'''
    (OUT/'notes.txt').write_text(notes,encoding='utf-8-sig')
    (OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>足関節・股関節 加速度波形</title><style>body{font:16px sans-serif;max-width:1300px;margin:30px auto;padding:16px;color:#183248}img{width:100%}pre{white-space:pre-wrap;background:#f3f6f9;padding:20px}h2{margin-top:48px}</style><h1>足関節・股関節由来の加速度波形</h1><pre>'+html.escape(notes)+'</pre>'+''.join(cards),encoding='utf-8')
    with zipfile.ZipFile(OUT.with_suffix('.zip'),'w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.iterdir(): z.write(p,arcname=OUT.name+'/'+p.name)
        z.write(Path(__file__),arcname=OUT.name+'/'+Path(__file__).name)
        z.write(Path(__file__).with_name('head_imu_analysis.py'),arcname=OUT.name+'/head_imu_analysis.py')
    print('Verified 13 recordings / 78000 samples: unchanged decimal strings; no additional processing.')
    for axis in ['AP','ML']:
        vals=[r[3] for r in checks if r[1]==axis]
        print(axis,'sum vs head RMSE range:',min(vals),max(vals),'mm/s^2')


if __name__=='__main__': main()
