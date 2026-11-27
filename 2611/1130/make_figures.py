from pathlib import Path
import csv, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'analysis_three_cycles_0p5'
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Yu Gothic','font.size':10,'axes.unicode_minus':False})
with (SOURCE/'frequency_responses.csv').open(encoding='utf-8-sig') as f:data=list(csv.DictReader(f))
fig,axes=plt.subplots(3,2,figsize=(10,9),layout='constrained',sharex=True)
for r,(s,n) in enumerate([('seguchi','瀬口'),('otsubo','大坪'),('takato','たかと')]):
    for kind,label,col in [('head','頭部IMU','#0072B2'),('com','COM推定','#D55E00')]:
        b=[x for x in data if x['subject']==s and x['output_source']==kind and x['input_order']=='0' and x['output_order']=='0']
        f=np.array([float(x['frequency_hz']) for x in b]);g=np.array([float(x['gain_db']) for x in b]);p=np.array([float(x['phase_deg']) for x in b])
        axes[r,0].semilogx(f,g,'o-',ms=4,color=col,label=label)
        # Break the visual guide at the wrapped-phase discontinuity.
        xx=[];yy=[]
        for k,(hz,v) in enumerate(zip(f,p)):
            if k and abs(v-p[k-1])>180:xx.append(np.nan);yy.append(np.nan)
            xx.append(hz);yy.append(v)
        axes[r,1].semilogx(xx,yy,'o-',ms=4,color=col,label=label)
    axes[r,0].set_title(n+'  ゲイン');axes[r,0].set_ylabel('変位比 [dB]')
    axes[r,1].set_title(n+'  位相');axes[r,1].set_ylabel('位相 [deg]');axes[r,1].set_ylim(-195,195);axes[r,1].set_yticks([-180,-90,0,90,180])
    for ax in axes[r]:ax.grid(which='both',alpha=.25)
axes[0,0].legend(ncol=2,fontsize=9)
for ax in axes[-1]:ax.set_xlabel('周波数 [Hz]')
fig.savefig(OUT/'gain_phase_displacement.png',dpi=200);plt.close(fig)
equations=[
 r'$R_x=F_y,\qquad N_y=M_x-h_0F_y,\qquad h_0=0.0123\ \mathrm{m}$',
 r'$x_{\mathrm{COM}}=\frac{N_y-cR_x}{m_b g},\qquad c=\frac{J_b}{m_b l_b}+L_f+l_b$',
 r'$m_b=0.978M,\quad J_b=0.0473MH^2,\quad l_b=0.513H,\quad L_f=0.039H$',
 r'$a_{\mathrm{head}}=\mathbf{e}_{\mathrm{AP}}^{\mathsf{T}}\left(R(t)\mathbf{a}_{\mathrm{sensor}}-g\mathbf{e}_{\mathrm{up}}\right)$',
 r'$y(t)=b_0+b_1(t-\bar{t})+\sum_{k=1}^{10}\left[a_k\cos(2\pi f_k t)+b_k\sin(2\pi f_k t)\right]+\varepsilon(t)$',
 r'$Y(f_k)=a_k-i b_k,\qquad s=i2\pi f,\qquad V=sX,\quad A=s^2X$',
 r'$H(f)=Y(f)/U(f),\quad G_{\mathrm{dB}}=20\log_{10}|H|,\quad \phi=\arg H$',
 r'$\Delta\phi=\mathrm{wrap}_{[-180^\circ,180^\circ)}\left(\phi_{\mathrm{COM}}-\phi_{\mathrm{head}}\right)$',
 r'$\Delta\phi_{\mathrm{model}}=-360f\tau,\qquad E(\tau)=\sqrt{\frac{1}{K}\sum_{k=1}^{K}\left[\mathrm{wrap}\left(\Delta\phi_k+360f_k\tau\right)\right]^2}$'
]
for k,eq in enumerate(equations):
    fig=plt.figure(figsize=(10,.62 if k!=1 else .8));fig.text(.01,.50,eq,fontsize=17,va='center',color='black')
    fig.savefig(OUT/f'equation_{k+1}.png',dpi=180,bbox_inches='tight',pad_inches=.07);plt.close(fig)
for name in ['seguchi','otsubo','takato']:
    shutil.copy2(SOURCE/f'{name}_head_com_bode_3x3.png',OUT/f'{name}_head_com_bode_3x3.png')
    shutil.copy2(SOURCE/f'{name}_phase_only.png',OUT/f'{name}_phase_only.png')
shutil.copy2(SOURCE/'head_com_constant_delay_fit.png',OUT/'constant_delay_fit.png')
for name in ['frequency_responses.csv','phase_difference_com_minus_head.csv','head_com_delay_fits.csv','method.json','head_com_delay_method.json']:
    dest=ROOT/'data';dest.mkdir(exist_ok=True);shutil.copy2(SOURCE/name,dest/name)
print('Figures and source values prepared')
