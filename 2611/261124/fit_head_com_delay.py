"""Exploratory circular phase fit; never interpret wrapped phase as unique lag."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

OUT=Path(__file__).resolve().parent/'three_trial_bode_sign_checked'
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
if args.output:OUT=args.output.resolve()
d=pd.read_csv(OUT/'phase_difference_com_minus_head.csv')
plt.rcParams.update({'font.family':'Yu Gothic','font.size':10,'axes.unicode_minus':False})
def wrap(x):return (x+180)%360-180
grid=np.linspace(-2,2,16001)
rows=[]
fig,axes=plt.subplots(3,2,figsize=(14,12),layout='constrained')
for row,(subject,name) in enumerate([('seguchi','瀬口'),('otsubo','大坪'),('takato','たかと')]):
    block=d[d.subject==subject].sort_values('frequency_hz')
    f=block.frequency_hz.to_numpy();phase=block.com_minus_head_deg.to_numpy()
    ax,profile=axes[row]
    ax.scatter(f,phase,c='#222222',s=38,zorder=4,label='実測：COM − 頭部')
    for offset,color,label in [(0,'#0072B2','純粋な遅れ'),(180,'#D55E00','符号反転＋遅れ')]:
        for band in ['all','low_0p7']:
            mask=np.ones(len(f),dtype=bool) if band=='all' else f<=.7000001
            ff,pp=f[mask],phase[mask]
            scores=np.mean(wrap(pp[:,None]-offset+360*ff[:,None]*grid)**2,axis=0)
            k=np.argmin(scores)
            opt=minimize_scalar(lambda t:np.mean(wrap(pp-offset+360*ff*t)**2),
                bounds=(grid[max(0,k-1)],grid[min(len(grid)-1,k+1)]),method='bounded')
            tau=float(opt.x);rmse=float(np.sqrt(opt.fun))
            rows.append(dict(subject=subject,band=band,offset_deg=offset,delay_s=tau,circular_rmse_deg=rmse,
                zero_delay_rmse_deg=float(np.sqrt(np.mean(wrap(pp-offset)**2))),search_min_s=-2,search_max_s=2,n_points=len(ff)))
            if band!='all':continue
            dense=np.linspace(float(f.min()),float(f.max()),1200);curve=wrap(offset-360*dense*tau)
            curve[np.r_[False,abs(np.diff(curve))>180]]=np.nan
            ax.plot(dense,curve,color=color,label=f'{label}：τ={tau:+.3f} s, RMSE={rmse:.1f}°')
            profile.plot(grid,np.sqrt(scores),color=color,label=label)
            profile.scatter([tau],[rmse],color=color,zorder=4)
    ax.axhline(0,c='gray',lw=.7);ax.set(title=name,xlabel='周波数 [Hz]',ylabel='位相差 [deg]',ylim=(-195,195),yticks=[-180,-90,0,90,180])
    ax.legend(fontsize=8,loc='best');ax.grid(alpha=.25)
    profile.axvline(0,c='gray',ls='--',lw=.8)
    profile.set(title=f'{name}：時間遅れを変えたときの当てはまり',xlabel='τ [s]：正＝COMが遅れる／負＝COMが先行',ylabel='円周上の位相誤差 RMSE [deg]',xlim=(-2,2))
    profile.legend(fontsize=8);profile.grid(alpha=.25)
fig.suptitle(f'頭部に対するCOMの位相差は、一定時間遅れで説明できるか？\n純粋な遅れ：Δφ=−360fτ ／ 符号反転＋遅れ：Δφ=180−360fτ（±180°表示）\n{len(f)}周波数を等重みで円周上に適合・探索範囲 ±2秒（探索上の仮定）・機器間遅延未校正',fontsize=13)
fig.savefig(OUT/'head_com_constant_delay_fit.png',dpi=150);plt.close(fig)
result=pd.DataFrame(rows);result.to_csv(OUT/'head_com_delay_fits.csv',index=False)
(OUT/'head_com_delay_method.json').write_text(json.dumps(dict(
    objective='Equal-weight mean squared principal circular phase residual; independent fits per participant.',
    models=['delta=-360*f*tau','delta=180-360*f*tau'],search_seconds=[-2,2],
    bands={'all':[float(f.min()),float(f.max())],'low_0p7':[float(f.min()),min(.7,float(f.max()))]},
    limitations=['Exploratory descriptive fit; no uncertainty or coherence available.',
    'Search range is an assumption, not a measured bound. Frequency grid 0.05 Hz gives exact delay aliasing every 20 seconds.',
    'Negative fitted delay means phase slope consistent with COM leading within this model, not causality.',
    'Estimated device timestamps, head gravity correction, COM axes/signs and anthropometry remain provisional.',
    'No free phase intercept; compare only zero and 180-degree intercepts. No manual phase unwrapping.']
),ensure_ascii=False,indent=2),encoding='utf-8')
print(result.round(3).to_string(index=False))
