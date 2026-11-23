"""Export head acceleration traces, one-sided FFT amplitude and RMS locally."""
from pathlib import Path
import csv
import io
import html
import zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('C:/Users/atvma/OneDrive/ドキュメント/MedicalForcePlate/csv')
NAMES = ['gridtrans05_seguchi00002', 'gridhip02_seguchi00003',
         'gridhip02_seguchi00004', 'gridhip05_seguchi00001',
         'gridhip05_seguchi00002', 'gridleg_seguchi00001',
         'gridleg02_seguchi00002', 'gridleg05_seguchi00001',
         'gridleg05_seguchi00002', 'gridstatics_seguchi00001',
         'gridtrans02_seguchi00001', 'gridtrans02_seguchi00002',
         'gridtrans05_seguchi00001']
OUT = Path(__file__).resolve().parent / 'outputs' / 'head_imu_20261026'


def fft_amplitude(x, fs):
    x = np.asarray(x) - np.mean(x)
    w = np.hanning(len(x))
    amplitude = np.abs(np.fft.rfft(x*w)) / w.sum()
    amplitude[1:-1 if len(x) % 2 == 0 else None] *= 2
    return np.fft.rfftfreq(len(x), 1/fs), amplitude


def write_csv(path, header, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def main(kind='head'):
    global OUT
    configs = {
        'head': ('頭部加速度', ['axhd', 'ayhd'], ['RMSHAap', 'RMSHAml']),
        'ankle': ('足関節モード加速度', ['xix1', 'xiy1'], ['RMSXI1ap', 'RMSXI1ml']),
        'hip': ('股関節モード加速度', ['xix2', 'xiy2'], ['RMSXI2ap', 'RMSXI2ml']),
    }
    title, columns, exported_columns = configs[kind]
    OUT = Path(__file__).resolve().parent / 'outputs' / f'{kind}_imu_20261026'
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams['font.family'] = ['Yu Gothic', 'Meiryo', 'DejaVu Sans']
    summary, cards, warnings = [], [], []
    for name in sorted(NAMES):
        path = BASE / f'261022_{name}.csv'
        raw = path.read_bytes()
        try:
            text = raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            text = raw.decode('cp932')
        rows = list(csv.reader(io.StringIO(text)))
        k = [i for i, row in enumerate(rows) if all(c in row for c in columns)][-1]
        header = rows[k]
        units = rows[[i for i, r in enumerate(rows[:k]) if len(r) == len(header) and r[29] == 'mm/s/s'][-1]]
        cols = [header.index(c) for c in columns]
        assert all(units[c] == 'mm/s/s' for c in cols)
        a = np.array([r for r in rows[k+1:] if r and any(v.strip() for v in r)], dtype=float)
        t, signals = a[:, 0], a[:, cols]
        dt = np.diff(t)
        assert np.isfinite(a).all() and (dt > 0).all()
        assert np.allclose(dt, np.median(dt), rtol=.001, atol=1e-8)
        fs = 1 / np.median(dt)
        assert np.isclose(fs, 100)
        mask = (t >= 5-1e-9) & (t < 35-1e-9)
        assert mask.sum() == 3000
        selected = signals[mask]
        centered = selected-selected.mean(axis=0)
        rms = np.sqrt(np.mean(centered**2, axis=0))
        full_rms = np.std(signals, axis=0)
        raw_rms = np.sqrt(np.mean(selected**2, axis=0))
        f, amp_ap = fft_amplitude(selected[:, 0], fs)
        _, amp_ml = fft_amplitude(selected[:, 1], fs)
        # Check FFT normalization independently using a bin-centered sinusoid.
        test_t = np.arange(3000)/fs
        tf, ta = fft_amplitude(3*np.sin(2*np.pi*test_t), fs)
        assert abs(ta[np.argmin(abs(tf-1))]-3) < 1e-4
        rect = np.fft.rfft(centered, axis=0)
        power = np.abs(rect)**2
        power[1:-1] *= 2
        assert np.allclose(power.sum(axis=0)/len(centered)**2, rms**2)
        export_header = next(r for r in rows[:k] if exported_columns[0] in r)
        export_row = rows[rows.index(export_header)+2]
        exported = [float(export_row[export_header.index(c)]) for c in exported_columns]
        if not np.allclose(rms, exported, rtol=1e-6, atol=1e-6):
            warnings.append(f'{path.name}: 装置内{exported_columns[0][:-2]}と5–35秒の平均除去RMSは一致しません。装置内指標を置換せず併記。')
        if rows[0][1] != '261022_'+name[:-5]:
            warnings.append(f'{path.name}: 内部測定名={rows[0][1]}。条件名はファイル名で表示。')
        for j, axis in enumerate(['AP', 'ML']):
            amp = [amp_ap, amp_ml][j]
            peak = 1 + np.argmax(amp[1:])
            summary.append([path.name, rows[0][1], axis, columns[j],
                            5,35,int(mask.sum()),fs,rms[j],raw_rms[j],selected[:,j].mean(),
                            full_rms[j],exported[j],f[peak],amp[peak]])
        stem = path.stem
        write_csv(OUT/f'{stem}_waveform.csv',
                  ['time_s',f'AP_{columns[0]}_mm_s2',f'ML_{columns[1]}_mm_s2','in_5_35s'],
                  zip(t,signals[:,0],signals[:,1],mask.astype(int)))
        write_csv(OUT/f'{stem}_fft.csv',
                  ['frequency_Hz','AP_amplitude_mm_s2','ML_amplitude_mm_s2'],zip(f,amp_ap,amp_ml))
        fig, axes = plt.subplots(2,2,figsize=(13,7),layout='constrained')
        fig.suptitle(f'{stem}\n{title} | 100 Hz | FFT・RMS: 5 ≤ t < 35 s',fontsize=14)
        for j,(label,color,amp) in enumerate([('AP 前後','#167d9a',amp_ap),('ML 左右','#d37828',amp_ml)]):
            ax=axes[j,0]
            ax.plot(t,signals[:,j],color=color,lw=.7)
            ax.axvspan(5,35,color=color,alpha=.1,label='FFT / RMS interval')
            ax.set(title=f'{label} 時刻波形（装置出力値）',xlabel='Time (s)',ylabel='Acceleration (mm/s²)',xlim=(0,60))
            ax.legend(fontsize=8)
            ax=axes[j,1]
            ax.plot(f,amp,color=color,lw=1)
            ax.set(title=f'{label} FFT | RMS = {rms[j]:.3f} mm/s²',xlabel='Frequency (Hz)',ylabel='One-sided amplitude (mm/s²)',xlim=(0,5))
            ax.text(.98,.95,f'全区間RMS: {full_rms[j]:.3f} mm/s²',ha='right',va='top',transform=ax.transAxes,fontsize=9)
        for ax in axes.flat:
            ax.grid(alpha=.2)
        fig.savefig(OUT/f'{stem}.png',dpi=150)
        plt.close(fig)
        cards.append(f'<h2>{html.escape(stem)}</h2><img src="{stem}.png"><p><a href="{stem}_waveform.csv">時刻波形CSV</a> · <a href="{stem}_fft.csv">FFT CSV（0–50 Hz）</a></p>')
    write_csv(OUT/'rms_summary.csv', ['file','internal_recording_name','axis','column','start_s','end_s_exclusive','samples','sampling_hz','RMS_demeaned_mm_s2','RMS_uncentered_mm_s2','mean_mm_s2','full_0_60s_RMS_demeaned_mm_s2','device_exported_RMS_mm_s2','FFT_peak_nonDC_Hz','FFT_peak_amplitude_mm_s2'],summary)
    notes = f'''{title}解析（13記録）
使用列: {columns[0]}=前後(AP)、{columns[1]}=左右(ML)、単位 mm/s²。装置処理済みの加速度。足関節・股関節モードは装置の推定値であり、関節角度や関節そのものの実測加速度ではありません。
全記録100 Hz、6000点、0–59.99秒。主解析は既存プロジェクトに合わせ5 ≤ t < 35秒（3000点）。
時刻波形は全区間の装置出力値、色付き背景は主解析区間。波形CSVも装置出力値。
主RMS=sqrt(mean((x-mean(x))²))。平均除去なしRMS、区間平均、全区間の平均除去RMS、装置内RMSもCSVに併記。
FFT: 区間平均除去→対称Hann窓→rFFT。片側振幅=abs(rFFT)/sum(window)、DC/Nyquist以外2倍。単位mm/s²。
周波数分解能1/30≈0.03333 Hz、Nyquist=50 Hz。図は0–5 Hz、FFT CSVは0–50 Hz。
RMSは窓掛け前の信号で算出。追加フィルタ、補間、再サンプリングなし。機器メタデータのフィルタ設定はHP=0.1 Hz、LP=4 Hz。
機器内の姿勢・重力補正は再実装していません。rawAx1/2等からのセンサー同定は行っていません。
ファイル名から推測した条件を内部名に上書きしません。元CSVは変更していません。
検証: 有限値・時刻連続性・サンプリング・単位確認、既知正弦波のFFT振幅確認、窓なしFFTのParseval等式によるRMS検証。
'''
    notes += '\n確認事項:\n'+'\n'.join(warnings)
    (OUT/'notes.txt').write_text(notes,encoding='utf-8-sig')
    body = '<!doctype html><meta charset="utf-8"><title>'+html.escape(title)+'解析</title><style>body{font:16px sans-serif;max-width:1200px;margin:32px auto;padding:16px;color:#173047}img{width:100%}pre{white-space:pre-wrap;background:#f2f5f8;padding:20px}h2{margin-top:48px}</style><h1>'+html.escape(title)+'：時刻波形・FFT・RMS</h1><p><a href="rms_summary.csv">RMS一覧CSV</a></p><pre>'+html.escape(notes)+'</pre>'+''.join(cards)
    (OUT/'index.html').write_text(body,encoding='utf-8')
    fig,ax=plt.subplots(figsize=(13,7),layout='constrained')
    labels=[f'261022_{n}' for n in sorted(NAMES)]
    ap=[r[8] for r in summary if r[2]=='AP']; ml=[r[8] for r in summary if r[2]=='ML']
    y=np.arange(len(labels))
    ax.barh(y-.18,ap,.36,label='AP 前後',color='#167d9a')
    ax.barh(y+.18,ml,.36,label='ML 左右',color='#d37828')
    ax.set(yticks=y,yticklabels=[s.replace('261022_','') for s in labels],xlabel='平均除去後 RMS (mm/s²)',title=f'{title} RMS | 5 ≤ t < 35 s')
    ax.invert_yaxis(); ax.legend(); ax.grid(axis='x',alpha=.2)
    fig.savefig(OUT/'rms_comparison.png',dpi=150); plt.close(fig)
    archive=OUT.with_suffix('.zip')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in OUT.iterdir():
            z.write(p,arcname=OUT.name+'/'+p.name)
        z.write(Path(__file__),arcname=OUT.name+'/head_imu_analysis.py')
    print(f'Completed: {len(NAMES)} files; {len(summary)} RMS rows; {OUT}')
    print('\n'.join(warnings).encode('ascii', errors='backslashreplace').decode('ascii'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['head', 'ankle', 'hip'], default='head')
    main(parser.parse_args().mode)
