from pathlib import Path
import csv
import math
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'outputs'/'joint_acceleration_20261101'
OUT=ROOT/'output'/'pdf'
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':['Yu Gothic','Meiryo','DejaVu Sans'],
                     'pdf.fonttype':42,'font.size':8})
files=sorted(SOURCE.glob('*_joint_acceleration.csv'))
assert len(files)==13
records=[]
for file in files:
    with file.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.reader(f))
    a=np.array(rows[1:],float)
    assert a.shape==(6000,5) and np.isfinite(a).all()
    a=a[a[:,0]>=5]
    assert a.shape==(5500,5) and a[0,0]==5
    label=file.name.removesuffix('_joint_acceleration.csv').removeprefix('261022_')
    if label=='gridleg_seguchi00001': label='gridleg02_seguchi00001'
    records.append((label,a))
records.sort(key=lambda record: record[0])
common=math.ceil(max(np.max(abs(a[:,1:3])) for _,a in records)/100)*100
colors=['#147b9a','#d36a22']
for mode in ['individual','common']:
    fig,axes=plt.subplots(4,4,figsize=(16.535,11.693))
    fig.subplots_adjust(left=.05,right=.985,bottom=.09,top=.875,wspace=.25,hspace=.55)
    fig.suptitle('足関節戦略・股関節戦略による頭部加速度 | AP方向・全13記録',
                 x=.05,y=.97,ha='left',fontsize=20,fontweight='bold',color='#193448')
    desc='記録ごとの縦軸：波形の形を比較（記録間で振幅の見た目は比較できません）' if mode=='individual' else f'共通の縦軸 ±{common:g} mm/s²：記録間の加速度の大きさを比較'
    fig.text(.05,.925,desc,fontsize=11)
    fig.legend(handles=[Line2D([0],[0],color=c,lw=1.8,label=l) for c,l in zip(colors,['足関節戦略 (xix1)','股関節戦略 (xix2)'])],
               loc='upper right',bbox_to_anchor=(.985,.94),ncol=2,frameon=False,fontsize=10)
    for ax,(label,a) in zip(axes.flat,records):
        for j,c in enumerate(colors):
            ax.plot(a[:,0],a[:,j+1],color=c,lw=.5,alpha=.85)
        limit=common if mode=='common' else max(25,math.ceil(np.max(abs(a[:,1:3]))/25)*25)
        ax.set(xlim=(5,60),ylim=(-limit,limit),xticks=[5,20,40,60],yticks=[-limit,0,limit])
        ax.set_title(label,loc='left',fontsize=9,pad=6,fontweight='bold')
        ax.tick_params(labelsize=7)
        ax.axhline(0,color='#89949b',lw=.4,zorder=0)
        ax.grid(alpha=.16,lw=.5)
        ax.spines[['top','right']].set_visible(False)
    for ax in list(axes.flat)[13:]: ax.set_axis_off()
    note_ax=list(axes.flat)[13]
    note_ax.text(0,1,'表示条件\n\n横軸：計測開始からの時刻 (s)\n縦軸：加速度 (mm/s²)\n最初の5秒を除外：5-59.99秒、100 Hz\n追加フィルタ・平均除去・正規化なし\ngridleg_seguchi00001 → gridleg02_seguchi00001',
                 transform=note_ax.transAxes,va='top',fontsize=9,linespacing=1.7)
    fig.text(.05,.045,'Sonobe & Mitsutake (2025), §2.2 | 装置算出済みのモード列を表示。生データからの独立再推定ではありません。',fontsize=8,color='#4c5b65')
    fig.text(.05,.022,'元CSVの装置側処理は残ります。gridtrans05の2件は内部測定名がgridtrans02です。条件表示はファイル名に従います。',fontsize=8,color='#4c5b65')
    target=OUT/f'AP_13recordings_after5s_{mode}_scale.pdf'
    fig.savefig(target)
    plt.close(fig)
    reader=PdfReader(target)
    assert len(reader.pages)==1
    text=reader.pages[0].extract_text()
    assert all(label in text for label,_ in records)
    assert 'gridleg_seguchi00001' not in [label for label,_ in records]
    desktop=Path.home()/'OneDrive'/'Pictures'/'Desktop'/'joint_acceleration_20261101'
    if desktop.is_dir(): shutil.copy2(target,desktop/target.name)
    print(target)
