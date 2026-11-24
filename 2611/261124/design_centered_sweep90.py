"""Centered 90s sweep with exact analytic acceleration envelope and rest velocities.

Solve for a slowly varying displacement envelope, rather than integrate from rest
and retain a secular position term. Slightly adjust the frequency schedule so
both endpoint velocities vanish. No onset, offset, or return interval is added.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp,quad,cumulative_simpson
from scipy.optimize import brentq
from design_visual_stimuli import A,metrics

OUT=Path(__file__).resolve().parent/'stimulus_design'/'sweep90_centered'
T,F0,F1=90.,.70,.05
L=np.log(F1/F0)

def angular(t,b):
    u=np.asarray(t)/T
    return 2*np.pi*F0*np.exp(L*u+b*u*(1-u))

def solve_envelope(b,tol=1e-8):
    def ode(t,y):
        w=angular(t,b);lam=(L+b*(1-2*t/T))/T;dlam=-2*b/T**2
        q=(2*y[1]-3*lam*y[0])/w
        return np.vstack([y[1],4*lam*y[1]-(4*lam**2-2*dlam)*y[0]+w*w*(y[0]-np.sqrt(np.maximum(1-q*q,1e-12)))])
    grid=np.linspace(0,T,901)
    sol=solve_bvp(ode,lambda a,z:np.array([a[1],z[1]]),grid,np.vstack([np.ones(len(grid)),np.zeros(len(grid))]),tol=tol,max_nodes=20000)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def end_phase(b):
    # y'=0 at both boundaries implies r'/r=-2*lambda there.
    lam0=(L+b)/T;lam1=(L-b)/T
    phi0=np.arctan2(-2*lam0,angular(0,b))
    phi1=np.arctan2(-2*lam1,angular(T,b))
    return quad(lambda t:float(angular(t,b)),0,T,epsabs=1e-11)[0]+phi0-phi1-44*np.pi

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    b=brentq(end_phase,-.5,.5,xtol=1e-13)
    sol=solve_envelope(b)
    t=np.arange(90001)/1000;w=angular(t,b);lam=(L+b*(1-2*t/T))/T;dlam=-2*b/T**2
    y,dy=sol.sol(t);r=A*y/w**2;dr=A*(dy-2*lam*y)/w**2
    q=(2*dy-3*lam*y)/w
    assert max(abs(q))<1
    # Exact displacement derivatives give acceleration magnitude A.
    ddr=r*w*w-A*np.sqrt(1-q*q)
    phi0=np.arctan2(dr[0],r[0]*w[0])
    phi=cumulative_simpson(w,x=t,initial=0)+phi0
    center=r[0]*np.cos(phi0)
    x=center-r*np.cos(phi)
    v=-dr*np.cos(phi)+r*w*np.sin(phi)
    coscoef=r*w*w-ddr;sincoef=2*dr*w+r*lam*w
    a=coscoef*np.cos(phi)+sincoef*np.sin(phi)
    acc_envelope=np.hypot(coscoef,sincoef)
    np.testing.assert_allclose(acc_envelope,A,rtol=1e-10,atol=1e-12)
    assert max(abs(v[[0,-1]]))<1e-8
    assert np.all(np.diff(w)<0) and np.isclose(w[-1]/(2*np.pi),F1)
    # Independent finite-difference check of exported position/velocity/acceleration.
    mid=slice(2,-2)
    assert max(abs(np.gradient(x,.001,edge_order=2)[mid]-v[mid]))<2e-6
    assert max(abs(np.gradient(v,.001,edge_order=2)[mid]-a[mid]))<2e-6
    upper=center+r;lower=center-r
    np.testing.assert_allclose((upper+lower)/2,center,atol=1e-14)
    freq=w/(2*np.pi)
    # Actual acceleration phase differs slightly from displacement phase in a chirp.
    acc_phase=phi-np.arctan2(sincoef,coscoef)
    acc_freq=np.gradient(acc_phase,t)/(2*np.pi)
    pd.DataFrame(dict(time_s=t,ap_displacement_m=x,ap_velocity_m_s=v,ap_acceleration_m_s2=a,displacement_carrier_frequency_hz=freq,acceleration_instantaneous_frequency_hz=acc_freq,displacement_center_m=np.full(len(t),center),displacement_upper_envelope_m=upper,displacement_lower_envelope_m=lower,acceleration_envelope_m_s2=acc_envelope)).to_csv(OUT/'centered_sweep90.csv',index=False)
    cfg=dict(status='design_only_not_installed_in_VR',duration_s=T,start_frequency_hz=F0,end_frequency_hz=F1,frequency_definition='displacement carrier instantaneous frequency; acceleration frequency differs slightly due to exact chirp envelope compensation',frequency_shape_b=float(b),frequency_formula='f(t)=0.70*exp(log(0.05/0.70)*u+b*u*(1-u)); u=t/90',acceleration_amplitude_m_s2=A,displacement_center_m=float(center),position_formula='x(t)=center-r(t)*cos(phi(t)); phi_dot=2*pi*f; r from BVP ensuring hypot(r*w^2-r_ddot,2*r_dot*w+r*w_dot)=A',envelope_equation='r=A*y/w^2; y_ddot=4*lambda*y_dot-(4*lambda^2-2*lambda_dot)*y+w^2*(y-sqrt(1-((2*y_dot-3*lambda*y)/w)^2)); y_dot(0)=y_dot(90)=0',boundary_phase_rule='integral(w dt)+atan2(-2*lambda0,w0)-atan2(-2*lambda90,w90)=44*pi',bvp_tolerance=1e-8,phase_rad=float(phi0),maxima=metrics(x,v,a),initial_state=[float(x[0]),float(v[0]),float(a[0])],final_state=[float(x[-1]),float(v[-1]),float(a[-1])],max_acceleration_envelope_error=float(max(abs(acc_envelope-A))),acceleration_frequency_at_endpoints_hz=[float(acc_freq[0]),float(acc_freq[-1])],extra_phases_s=0,after_trial='Hold final position without reset; no return interval',notes='Fixed displacement envelope center, not zero moving-average over an arbitrary window. x starts at zero and endpoints have zero velocity. Acceleration steps at trial boundaries. No drift subtraction. Frequency shape slightly adjusted from exponential. Carrier spans 0.70..0.05Hz; acceleration instantaneous frequency is explicitly exported.')
    (OUT/'centered_sweep90.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(4,1,figsize=(11,9),sharex=True)
    axes[0].plot(t[::10],x[::10],lw=1)
    axes[0].plot(t[::50],upper[::50],'--',c='gray',lw=.7,label='Symmetric displacement envelope')
    axes[0].plot(t[::50],lower[::50],'--',c='gray',lw=.7)
    axes[0].axhline(center,color='black',ls=':',lw=1,label='Fixed center');axes[0].legend(fontsize=8)
    axes[1].plot(t[::10],v[::10],lw=.9)
    axes[2].plot(t[::10],a[::10],lw=.9)
    axes[2].axhline(A,color='gray',ls=':',lw=.6);axes[2].axhline(-A,color='gray',ls=':',lw=.6)
    axes[3].plot(t[::10],freq[::10],label='Displacement carrier')
    axes[3].plot(t[::10],acc_freq[::10],ls='--',label='Acceleration instantaneous frequency');axes[3].legend(fontsize=8)
    for ax,label in zip(axes,['Displacement [m]','Velocity [m/s]','Acceleration [m/s2]','Frequency [Hz]']):ax.set(ylabel=label,xlim=(0,90));ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time [s]');fig.suptitle('Centered 90s downsweep: fixed position center, constant acceleration envelope, no extra phases')
    fig.tight_layout();fig.savefig(OUT/'centered_sweep90.png',dpi=170);plt.close(fig)
    print(json.dumps(cfg,indent=2))

if __name__=='__main__':main()
