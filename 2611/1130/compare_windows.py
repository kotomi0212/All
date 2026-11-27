from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
old=R.parent/'261124/three_trial_bode_sign_checked';new=R/'analysis_three_cycles_0p5'
a=pd.read_csv(old/'phase_difference_com_minus_head.csv');b=pd.read_csv(new/'phase_difference_com_minus_head.csv')
m=a.merge(b,on=['subject','frequency_hz'],suffixes=('_4cycles','_3cycles'))
m['phase_change_deg']=(m.com_minus_head_deg_3cycles-m.com_minus_head_deg_4cycles+180)%360-180
m.to_csv(R/'data/window_phase_comparison.csv',index=False)
plt.rcParams.update({'font.family':'Yu Gothic','font.size':10,'axes.unicode_minus':False})
fig,axes=plt.subplots(3,2,figsize=(12,10),layout='constrained')
g0=pd.read_csv(old/'frequency_responses.csv');g1=pd.read_csv(new/'frequency_responses.csv')
for k,(s,n) in enumerate([('seguchi','瀬口'),('otsubo','大坪'),('takato','たかと')]):
 block=m[m.subject==s]
 for key,label,col in [('4cycles','4周期 10–90秒','#999999'),('3cycles','3周期 30–90秒','#0072B2')]:
  axes[k,0].semilogx(block.frequency_hz,block[f'com_minus_head_deg_{key}'],'o',label=label,color=col)
 axes[k,0].set(title=n+'  COM−頭部の位相差',ylim=(-195,195),ylabel='位相差 [deg]');axes[k,0].axhline(0,color='gray',lw=.6)
 for kind,label,col in [('head','頭部IMU','#0072B2'),('com','COM推定','#D55E00')]:
  mask=lambda d:(d.subject==s)&(d.output_source==kind)&(d.input_order==0)&(d.output_order==0)
  x=g0[mask(g0)].merge(g1[mask(g1)],on='frequency_hz',suffixes=('_old','_new'))
  axes[k,1].semilogx(x.frequency_hz,x.gain_db_new-x.gain_db_old,'o-',label=label,color=col)
 axes[k,1].set(title=n+'  ゲインの変化 3周期−4周期',ylabel='ゲイン変化 [dB]');axes[k,1].axhline(0,color='gray',lw=.6)
 for ax in axes[k]:ax.grid(which='both',alpha=.25);ax.legend(fontsize=8);ax.set_xlabel('周波数 [Hz]')
 print(s,'median phase change',round(abs(block.phase_change_deg).median(),1),'points >30 deg',sum(abs(block.phase_change_deg)>30))
fig.suptitle('1周期目の除外による解析区間への依存性',fontsize=15)
fig.savefig(R/'figures/window_comparison.png',dpi=160)
