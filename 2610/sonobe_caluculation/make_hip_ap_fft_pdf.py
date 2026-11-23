"""AP hip FFT, 5 <= t < 60 s, rectangular window, display <= 2 Hz."""
from pathlib import Path
import csv
import math
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'outputs'/'joint_acceleration_20261101'
OUT=ROOT/'output'/'pdf'
DATA=ROOT/'outputs'/'hip_AP_fft_after5s'
OUT.mkdir(parents=True,exist_ok=True)
DATA.mkdir(parents=True,exist_ok=True)
DESKTOP=Path.home()/'OneDrive'/'Pictures'/'Desktop'/'joint_acceleration_20261101'
plt.rcParams.update({'font.family':['Yu Gothic','Meiryo','DejaVu Sans'],'pdf.fonttype':42,'font.size':8})
records=[]
summary=[]
for path in sorted(SOURCE.glob('*_joint_acceleration.csv')):
    with path.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        rows=list(reader)
    t=np.array([float(r['time_s']) for r in rows])
    x=np.array([float(r['AP_hip_xix2_mm_s2']) for r in rows])
    assert len(t)==6000 and np.allclose(np.diff(t),.01,atol=1e-8)
    x=x[(t>=5)&(t<60)]
    assert len(x)==5500 and np.isfinite(x).all()
    label=path.name.removesuffix('_joint_acceleration.csv').removeprefix('261022_')
    if label=='gridleg_seguchi00001': label='gridleg02_seguchi00001'
    centered=x-x.mean()
    fft=np.fft.rfft(centered)
    freq=np.fft.rfftfreq(len(x),1/100)
    amplitude=np.abs(fft)/len(x)
    amplitude[1:-1]*=2
    power=np.abs(fft)**2
    power[1:-1]*=2
    assert np.isclose(power.sum()/len(x)**2,np.mean(centered**2))
    # Rectangular-window amplitude must recover a bin-centered reference sinusoid.
    test=np.sin(2*np.pi*(55/55)*np.arange(len(x))/100)*3
    assert abs(2*abs(np.fft.rfft(test)[55])/len(x)-3)<1e-10
    mask=freq<=2+1e-12
    f,a=freq[mask],amplitude[mask]
    assert len(f)==111 and np.isclose(f[-1],2)
    peak=1+np.argmax(a[1:])
    summary.append([label,5,60,len(x),100,1/55,f[peak],a[peak]])
    records.append((label,f,a))
assert len(records)==13
records.sort(key=lambda r:('statics' in r[0],r[0]));summary.sort(key=lambda r:('statics' in r[0],r[0]))
csv_path=DATA/'AP_hip_FFT_0_2Hz_after5s.csv'
with csv_path.open('w',newline='',encoding='utf-8-sig') as stream:
    w=csv.writer(stream)
    w.writerow(['frequency_Hz']+[r[0]+'_amplitude_mm_s2' for r in records])
    w.writerows(zip(records[0][1],*[r[2] for r in records]))
with (DATA/'peak_summary.csv').open('w',newline='',encoding='utf-8-sig') as stream:
    w=csv.writer(stream); w.writerow(['recording','start_s','end_s_exclusive','samples','sampling_hz',
                                     'frequency_resolution_Hz','peak_frequency_0_2Hz_nonDC','peak_amplitude_mm_s2'])
    w.writerows(summary)
max_amp=max(float(r[2].max()) for r in records)
common=math.ceil(max_amp/2)*2
for mode in ['common']:
    fig,axes=plt.subplots(4,4,figsize=(16.535,11.693))
    fig.subplots_adjust(left=.05,right=.985,bottom=.09,top=.875,wspace=.25,hspace=.55)
    fig.suptitle('股関節戦略による頭部加速度のFFT | AP方向・0-2 Hz・全13記録',
                 x=.05,y=.97,ha='left',fontsize=19,fontweight='bold',color='#193448')
    desc='記録ごとの縦軸：周波数分布の形を比較' if mode=='individual' else '共通の縦軸：周波数成分の振幅を記録間で比較'
    fig.text(.05,.925,desc+' | 解析区間：5-59.99秒',fontsize=11)
    for ax,(label,f,a) in zip(axes.flat,records):
        ax.plot(f,a,color='#d36a22',lw=.9)
        marker_hz = .2 if any(s in label for s in ['hip02','leg02','trans02']) else .5 if any(s in label for s in ['hip05','leg05','trans05']) else None
        if marker_hz is not None:
            ax.axvline(marker_hz,color='#2464a0',ls='--',lw=1.1,zorder=3)
            ax.text(marker_hz+.025,.96,f'{marker_hz:g} Hz',transform=ax.get_xaxis_transform(),
                    va='top',ha='left',fontsize=7,color='#2464a0',
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1))
        limit=common if mode=='common' else max(2,math.ceil(float(a.max())/2)*2)
        ax.set(xlim=(0,2),ylim=(0,limit),xticks=[0,.5,1,1.5,2])
        ax.set_title(label,loc='left',fontsize=9,pad=6,fontweight='bold')
        ax.tick_params(labelsize=7)
        ax.grid(alpha=.2,lw=.5)
        ax.spines[['top','right']].set_visible(False)
    for ax in list(axes.flat)[13:]: ax.set_axis_off()
    list(axes.flat)[13].text(0,1,'FFT条件\n\n入力：AP股関節モード (xix2)\n横軸：周波数 (Hz)\n縦軸：片側振幅 (mm/s²)\n100 Hz / 5500点 / 分解能 0.01818 Hz\n区間平均を除去、矩形窓（窓掛けなし）\n追加フィルタ・平滑化・補間なし',
        transform=list(axes.flat)[13].transAxes,va='top',fontsize=9,linespacing=1.6)
    fig.text(.05,.065,'青破線：hip02・leg02・trans02は0.2 Hz、hip05・leg05・trans05は0.5 Hz（ファイル名による条件表示）。',fontsize=8,color='#2464a0')
    fig.text(.05,.045,'2 Hz以下は表示範囲です。低域通過フィルタはかけていません。片側振幅 = |FFT| / N、DC・Nyquist以外は2倍。',fontsize=8,color='#4c5b65')
    fig.text(.05,.022,'装置算出済みモード列を使用。装置側の既存処理は残ります。矩形窓のため、周波数がFFTビンに一致しない成分には漏れが生じます。',fontsize=8,color='#4c5b65')
    target=OUT/f'AP_hip_FFT_0_2Hz_after5s_{mode}_scale_statics_last_marked.pdf'
    fig.savefig(target);plt.close(fig)
    pdf=PdfReader(target)
    assert len(pdf.pages)==1 and all(r[0] in pdf.pages[0].extract_text() for r in records)
    if DESKTOP.is_dir():shutil.copy2(target,DESKTOP/target.name)
    print(target)
if DESKTOP.is_dir():
    shutil.copy2(csv_path,DESKTOP/csv_path.name)
print('Verified: 13 records, 5500 samples each; FFT amplitude and Parseval checks passed.')
