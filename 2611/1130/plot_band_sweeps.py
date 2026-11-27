from pathlib import Path
import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent;out=root/'sweep100'
plt.rcParams.update({'font.family':'Yu Gothic','font.size':10,'axes.unicode_minus':False})
overview,axs=plt.subplots(3,2,figsize=(12,9),layout='constrained')
for band in range(1,4):
 for col,direction in enumerate(['up','down']):
  name=f'study_sweep100_band{band}_{direction}'
  d=pd.read_csv(out/f'{name}.csv');c=json.loads((out/f'{name}.json').read_text())
  title=f"帯域{band} {'上昇' if direction=='up' else '下降'} {c['start_frequency_hz']:.3f} → {c['end_frequency_hz']:.3f} Hz"
  fig,axes=plt.subplots(4,1,figsize=(11,8),sharex=True,layout='constrained')
  for ax,key,label in zip(axes,['displacement_m','velocity_m_s','acceleration_m_s2','frequency_hz'],['変位 [m]','速度 [m/s]','加速度 [m/s²]','周波数 [Hz]']):
   ax.plot(d.time_s,d[key],lw=.9,color='#0072B2');ax.set_ylabel(label);ax.grid(alpha=.25);ax.set_xlim(0,100)
   ax.axvspan(0,10,color='#cccccc',alpha=.4);ax.axvspan(95,100,color='#cccccc',alpha=.4)
  axes[-1].set_xlabel('時間 [s]');fig.suptitle(title+'\n立ち上げ10秒・掃引85秒・停止5秒 | 灰色は解析除外',fontsize=13)
  fig.savefig(out/f'{name}.png',dpi=160);plt.close(fig)
  ax=axs[band-1,col];ax.plot(d.time_s,d.displacement_m,lw=1,color='#0072B2');ax.set_title(title,fontsize=11);ax.set_ylabel('変位 [m]');ax.set_xlabel('時間 [s]');ax.set_xlim(0,100);ax.grid(alpha=.25)
  ax.axvspan(0,10,color='#cccccc',alpha=.4);ax.axvspan(95,100,color='#cccccc',alpha=.4)
overview.suptitle('100秒sweep 6パターン | 立ち上げ10秒・本体85秒・停止5秒\n本体の速度振幅0.10 m/s・停止後は終了位置を保持',fontsize=14)
overview.savefig(root/'figures/sweep100_six_patterns.png',dpi=180)
overview.savefig(root/'figures/sweep100_six_patterns.pdf')
print('Saved six detailed plots and overview')
