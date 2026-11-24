"""Gain-only Bode-style figure. Run with the repository .venv Python."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

OUT=Path(__file__).resolve().parent/'grid_com_identification'
plt.rcParams.update({'font.family':'Yu Gothic','font.size':11,'axes.unicode_minus':False})


def load(name):
    with (OUT/name).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))


def q(f):
    H=1.48; M=43.4; mb=.978*M; jb=.0481*M*H**2; lb=.504*H
    inertia=jb/lb+mb*lb; kp=10.53*inertia; kd=3.97*inertia
    s=2j*np.pi*np.asarray(f); delay=np.exp(-.155*s)
    return abs(-kp*delay/((kp+kd*s)*delay-mb*9.80665))


def curve(f,row):
    status=row['status']
    if status=='T_infinite_limit':
        return q(f)*float(row['gamma_over_T'])/(2*np.pi*f)
    return q(f)*float(row['gamma_s'])/np.sqrt(1+(2*np.pi*f*float(row['T_s']))**2)


def main():
    data=[r for r in load('com_harmonic_amplitudes.csv') if r['output']=='eq3_Fx_My']
    fits=[r for r in load('com_parameter_fits.csv') if r['output']=='eq3_Fx_My' and r['selection']=='all_four_trials']
    mainfit=next(r for r in fits if r['start_s']=='5' and r['end_s']=='35')
    primary=[r for r in data if r['start_s']=='5' and r['end_s']=='35']
    f=np.geomspace(.1,1,500)
    db=lambda x:20*np.log10(x)  # Numerical gain referenced to 1 second.
    fig,axes=plt.subplots(1,2,figsize=(13,5.6),layout='constrained')
    ax=axes[0]
    ax.semilogx(f,db(curve(f,mainfit)),color='#d55e00',lw=2,
                label='モデルの極限曲線（T→∞、有限同定値なし）')
    for trial,marker,color in [('00001','o','#0072b2'),('00002','s','#009e73')]:
        block=[r for r in primary if r['trial'].endswith(trial)]
        ax.scatter([float(r['frequency_Hz']) for r in block],
                   [db(float(r['measured_gain_s'])) for r in block],s=65,
                   marker=marker,color=color,label=f'実測：試行{int(trial)}',zorder=4)
    ax.set_title('主解析：記録開始後 5–35秒')
    ax.legend(fontsize=9,loc='lower left')
    rows=[]
    for row in fits:
        label=f"{row['start_s']}–{row['end_s']}秒"
        if row['status']=='finite':
            label+=f"：T={float(row['T_s']):.3f} s"
            style='-'
        elif row['status']=='T_zero_limit': label+='：T=0 の境界'; style=':'
        else: label+='：T→∞ の極限'; style='--'
        line,=axes[1].semilogx(f,db(curve(f,row)),ls=style,lw=1.7,label=label)
        block=[r for r in data if r['start_s']==row['start_s'] and r['end_s']==row['end_s']]
        for freq in [.2,.5]:
            mean=np.mean([float(r['measured_gain_s']) for r in block if float(r['frequency_Hz'])==freq])
            axes[1].scatter(freq,db(mean),s=30,color=line.get_color(),zorder=3)
        for freq,gain in zip(f,curve(f,row)):
            rows.append([row['start_s'],row['end_s'],row['status'],freq,gain,db(gain)])
    axes[1].set_title('区間による変化：点は2試行の平均ゲイン')
    axes[1].legend(fontsize=9,loc='lower left')
    for ax in axes:
        ax.set(xlabel='周波数 [Hz]',ylabel='COM変位 / 視覚速度ゲイン [dB re 1 s]',xlim=(.1,1),ylim=(-90,-15))
        ax.set_xticks([.1,.2,.5,1],labels=['0.1','0.2','0.5','1.0'])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.grid(which='both',alpha=.25)
        for ff in [.2,.5]: ax.axvline(ff,color='gray',alpha=.2,lw=.8)
    fig.suptitle('COM応答のボード線図（ゲインのみ）',fontsize=16)
    fig.supxlabel('実測は0.2・0.5 Hzの2点のみ。曲線は既存モデルによる補間・外挿。刺激との位相は未同期のため表示しない。',fontsize=10)
    for ext in ['png','svg']:
        fig.savefig(OUT/f'com_bode_gain.{ext}',dpi=180)
    plt.close(fig)
    with (OUT/'com_bode_model_curves.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(['start_s','end_s','status','frequency_Hz','gain_s','gain_dB_re_1_s']);writer.writerows(rows)
    print(OUT/'com_bode_gain.png')


if __name__=='__main__': main()
