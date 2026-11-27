"""Phase-only comparison and circular COM-minus-head phase differences."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent/'three_trial_bode_sign_checked'
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
if args.output:OUT=args.output.resolve()
import json
window=json.loads((OUT/'method.json').read_text(encoding='utf-8'))['window_s']
d=pd.read_csv(OUT/'frequency_responses.csv')
plt.rcParams.update({'font.family':'Yu Gothic','font.size':10,'axes.unicode_minus':False})
orders=['変位','速度','加速度']
rows=[]
for subject,name in [('seguchi','瀬口'),('otsubo','大坪'),('takato','たかと')]:
    fig,axes=plt.subplots(3,3,figsize=(15,10),sharex=True,sharey=True,layout='constrained')
    reference=None
    for j in range(3):
        for i in range(3):
            ax=axes[j,i]; arrays={}
            for kind,label,color in [('head','頭部IMU','#0072B2'),('com','COM推定','#D55E00')]:
                b=d[(d.subject==subject)&(d.output_source==kind)&(d.input_order==i)&(d.output_order==j)].sort_values('frequency_hz')
                arrays[kind]=b.phase_deg.to_numpy()
                f=b.frequency_hz.to_numpy()
                # Markers show wrapped phases; connecting lines are visual guides.
                ax.semilogx(f,b.phase_deg,'o-',ms=4,label=label,color=color)
            delta=(arrays['com']-arrays['head']+180)%360-180
            if reference is None:
                reference=delta
                for hz,h,c,v in zip(f,arrays['head'],arrays['com'],delta):
                    rows.append(dict(subject=subject,frequency_hz=hz,head_phase_deg=h,com_phase_deg=c,com_minus_head_deg=v))
            else:
                assert np.allclose(np.angle(np.exp(1j*np.deg2rad(delta-reference))),0,atol=1e-9)
            ax.set_title(f'刺激{orders[i]} → {orders[j]}')
            ax.set_ylim(-195,195);ax.set_yticks([-180,-90,0,90,180]);ax.grid(which='both',alpha=.25)
            if i==0:ax.set_ylabel('位相 [deg]')
            if j==2:ax.set_xlabel('周波数 [Hz]')
    axes[0,0].legend(ncol=2,fontsize=9)
    fig.suptitle(f'{name} | 頭部IMUとCOM推定：位相のみ\n2026/11/24・{window[0]:g}–{window[1]:g} s | 位相は±180°表示・機器間遅延未校正',fontsize=14)
    fig.savefig(OUT/f'{subject}_phase_only.png',dpi=140);plt.close(fig)
table=pd.DataFrame(rows)
table.to_csv(OUT/'phase_difference_com_minus_head.csv',index=False)
print(table.pivot(index='frequency_hz',columns='subject',values='com_minus_head_deg').round(1).to_string())
