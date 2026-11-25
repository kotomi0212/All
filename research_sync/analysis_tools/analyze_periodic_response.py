"""Cycle-resolved VR standing response audit. Read-only inputs; private outputs."""
import argparse,csv,json,hashlib,math
from pathlib import Path
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TARGET='study_multicosine10_equal_velocity'
def write(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
def table(path):
    with path.open(encoding='utf-8-sig',newline='') as f:r=list(csv.DictReader(f))
    if not r:raise ValueError('Empty CSV')
    data={}
    for k in r[0]:
        try:data[k]=np.array([float(x[k]) if x[k] else np.nan for x in r])
        except ValueError:pass
    return data,list(r[0])
def harmonic(t,y,f):
    X=np.column_stack([np.ones(len(t)),t-np.mean(t),*[a for hz in f for a in (np.cos(2*np.pi*hz*t),np.sin(2*np.pi*hz*t))]])
    if len(t)<X.shape[1]*3 or not np.isfinite(y).all():raise ValueError('Insufficient/nonfinite data')
    b=np.linalg.lstsq(X,y,rcond=None)[0]
    return b[2::2]-1j*b[3::2]
def period(f):
    frac=[Fraction(float(x)).limit_denominator(100000) for x in f]
    den=math.lcm(*(x.denominator for x in frac));num=math.gcd(*(x.numerator*(den//x.denominator) for x in frac))
    return den/num
def bh(p):
    order=np.argsort(p);v=np.array(p)[order]*len(p)/np.arange(1,len(p)+1)
    v=np.minimum.accumulate(v[::-1])[::-1];result=np.empty(len(p));result[order]=np.minimum(v,1);return result
def metrics(U,Y):
    n=len(U);H=np.sum(np.conj(U)*Y,axis=0)/np.sum(abs(U)**2,axis=0)
    Z=Y/U;mean=Z.mean(axis=0)
    variance=np.sum(abs(Z-mean)**2,axis=0)/(n-1) if n>1 else np.full(U.shape[1],np.nan)
    snr=10*np.log10(np.maximum(abs(mean)**2,1e-300)/variance)
    adjusted=abs(mean)**2-variance/n
    snr_adj=np.where(adjusted>0,10*np.log10(np.maximum(adjusted,1e-300)/variance),np.nan)
    coherence=abs(np.sum(np.conj(U)*Y,axis=0))**2/(np.sum(abs(U)**2,axis=0)*np.sum(abs(Y)**2,axis=0)) if n>1 else np.full(U.shape[1],np.nan)
    coherence=np.clip(coherence,0,1)
    p=(1-coherence)**(n-1)
    resultant=abs(np.mean(Z/np.maximum(abs(Z),1e-300),axis=0))
    phase_sd=np.degrees(np.sqrt(-2*np.log(np.minimum(resultant,1))))
    return H,Z,snr,snr_adj,coherence,p,bh(p),phase_sd

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();root=args.root.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    # Both historical datasets and acquisition copies are inventoried; duplicate trial IDs are never pooled.
    paths=list((root/'research_sync/recordings').rglob('*metadata.json'))+list((root/'2611').rglob('*metadata.json'))
    records=[]
    for p in paths:
        try:
            m=json.loads(p.read_text(encoding='utf-8-sig'));stim=m.get('stimulus',{})
            if not isinstance(stim,dict):continue
            if not m.get('trial_id'):continue
            records.append((p,m,stim))
        except (ValueError,OSError):continue
    main_records=[r for r in records if '/saved/' in r[0].as_posix() and r[2].get('condition')==TARGET]
    main_records.sort(key=lambda r:r[1].get('pc_commands',[{}])[0].get('pc_epoch_s',0))
    subjects=list(dict.fromkeys(r[1]['participant']['subject'] for r in main_records))
    aliases={s:'Participant '+chr(65+i) for i,s in enumerate(subjects)}
    inventory=[];used=set();selected=[]
    # Prefer named saved files over raw acquisition copies.
    records.sort(key=lambda r:('/saved/' not in r[0].as_posix(),str(r[0])))
    for p,m,stim in records:
        tid=m['trial_id'];cond=stim.get('condition','unknown');sub=(m.get('participant') or {}).get('subject','unknown')
        duplicate=tid in used;used.add(tid)
        target=cond==TARGET and not duplicate and m.get('termination')=='complete'
        status='candidate' if target else 'excluded'
        reason='' if target else ('duplicate trial_id (acquisition copy)' if duplicate else ('incomplete target trial; no full steady interval' if cond==TARGET else 'different condition; not pooled'))
        inventory.append(dict(metadata=str(p.relative_to(root)),trial_id=tid,participant=aliases.get(sub,'other'),condition=cond,termination=m.get('termination'),status=status,reason=reason,analysis_start_s='',analysis_end_s='',cycles=0))
        if target:selected.append((p,m,stim,inventory[-1]))
    results=[];reps=[];schema=[];inputs=[];sync=[];provenance=[];failures=[];subject_summary=[]
    fig_input,ax_input=plt.subplots(max(1,len(selected)),1,figsize=(13,3*max(1,len(selected))),squeeze=False,layout='constrained')
    for index,(mp,m,c,inv) in enumerate(selected):
      try:
        folder=mp.parent;stem=mp.name.removesuffix('_metadata.json')
        qp=folder/(stem+'_quest.csv');fp=folder/(stem+'_combalance.csv')
        if not qp.exists():qp=folder/'quest.csv';fp=folder/'combalance.csv'
        q,qcols=table(qp);b,bcols=table(fp)
        alias=inv['participant'];tid=m['trial_id'];f=np.array(c['frequencies_hz']);T=period(f)
        if len(f)!=10:raise ValueError('Not ten frequencies')
        start,stop=c['analysis_interval_s'];cycles=int(np.floor((stop-start)/T+1e-9));stop=start+cycles*T
        if m.get('termination')!='complete':raise ValueError('Incomplete trial')
        qt,qe,ft=q['pc_mapped_time_s'],q['stimulus_elapsed_s'],b['pc_estimated_sample_time_s']
        for label,t in [('Quest PC',qt),('Quest elapsed',qe),('force PC',ft)]:
            if not np.isfinite(t).all() or np.any(np.diff(t)<=0):raise ValueError(label+' invalid time')
        bounds=np.interp([start,stop],qe,qt)
        if qe[0]>start or qe[-1]<stop or bounds[0]<ft[0] or bounds[1]>ft[-1]:raise ValueError('Missing full analysis coverage')
        if max(np.max(np.diff(qt)),np.max(np.diff(ft)))>=.1:raise ValueError('Gap >=100 ms')
        # No shift optimization: invert recorded Quest-to-PC timestamp mapping.
        elapsed=np.interp(ft,qt,qe,left=np.nan,right=np.nan)
        mass=m['weight_estimate']['body_mass_kg'];height=m['participant']['height_cm']/100
        if not mass or mass<=0:raise ValueError('Missing mass')
        outputs={'COM_acceleration_Fy_over_M':(b['sensor_Fy_N']/mass,'m/s^2','AP Fy sign provisional; mass static Fz estimate'),
                 'Head_Z_raw':(b['external_imu_count_2'],'counts','sensor Z, user AP; not lab acceleration; no gravity correction')}
        # Previously plotted male/female curves used unassigned population coefficients.
        # Do not reconstruct subject COM position without confirmed body parameters.
        for hz in f:
            results.append(dict(participant=alias,trial_id=tid,output='COM_position_Sonobe',frequency_hz=hz,
                                gain='',phase_deg='',cycles=cycles,independent_trials=1,coherence='',snr_db='',
                                status='unavailable',reason='No confirmed participant lb/Jb; prior male/female defaults are assumptions; moment origin unvalidated'))
        for origin,cols,data in [('Quest',qcols,q),('ComBalance',bcols,b)]:
            for col in cols:
                values=data.get(col);unit='unknown'
                if col.endswith('_s'):unit='s'
                if col.endswith('_N'):unit='N'
                if col.endswith('_Nm'):unit='N m'
                if col.startswith('external_imu_count_'):unit='raw counts'
                if col.startswith('onboard_imu_count_') or col.startswith('force_adc_'):unit='raw counts'
                if col=='body_mass_kg':unit='kg'
                if col=='height_cm':unit='cm'
                if col=='sampling_hz':unit='Hz'
                if col=='ap_displacement_m':unit='m'
                if col=='ap_velocity_m_s':unit='m/s'
                if col=='ap_acceleration_m_s2':unit='m/s^2'
                schema.append(dict(participant=alias,trial_id=tid,source=origin,column=col,unit=unit,missing_numeric=int(np.count_nonzero(~np.isfinite(values))) if values is not None else 'not_numeric_or_not_loaded'))
        U=[];Ys={k:[] for k in outputs}
        grid=np.arange(start,stop,.01)
        uniform_velocity=np.interp(grid,qe,q['ap_velocity_m_s'])
        fft_f=np.fft.rfftfreq(len(grid),.01);fft_a=2*abs(np.fft.rfft(uniform_velocity-uniform_velocity.mean()))/len(grid)
        peak_indices=np.argsort(fft_a[1:])[-10:]+1
        observed_frequencies=np.sort(fft_f[peak_indices])
        if not np.allclose(observed_frequencies,np.sort(f),atol=1e-6):raise ValueError('Configuration and strongest ten logged input frequencies differ')
        write(out/(alias.replace(' ','_')+'_input_frequencies.csv'),[dict(config_hz=a,logged_fft_peak_hz=b,source_used='metadata exact frequencies; confirmed by logged velocity') for a,b in zip(f,observed_frequencies)])
        for cycle in range(cycles):
            a=start+cycle*T;z=a+T;qm=(qe>=a)&(qe<z);bm=(elapsed>=a)&(elapsed<z)
            if qm.sum()<100 or bm.sum()<100:raise ValueError('Insufficient samples in cycle')
            # +AP is forward (-world Z): world.position.z = base - logged x, so logged velocity is BACKGROUND +AP.
            u=harmonic(qe[qm],q['ap_velocity_m_s'][qm],f);U.append(u)
            X=harmonic(qe[qm],q['ap_displacement_m'][qm],f);A=harmonic(qe[qm],q['ap_acceleration_m_s2'][qm],f)
            for j,hz in enumerate(f):inputs.append(dict(participant=alias,cycle=cycle+1,frequency_hz=hz,input_velocity_amplitude=abs(u[j]),configured_velocity_amplitude=c.get('velocity_component_m_s'),relative_error_vs_config=abs(abs(u[j])/c['velocity_component_m_s']-1),velocity_vs_jomegaX_relative_error=abs(u[j]-2j*np.pi*hz*X[j])/abs(u[j]),acceleration_vs_jomegaV_relative_error=abs(A[j]-2j*np.pi*hz*u[j])/abs(A[j]),quest_samples=int(qm.sum()),force_samples=int(bm.sum())))
            for k,(y,unit,note) in outputs.items():Ys[k].append(harmonic(elapsed[bm],y[bm],f))
        U=np.array(U)
        fig,axes=plt.subplots(4,len(outputs),figsize=(12,12),sharex=True)
        fig.subplots_adjust(left=.11,right=.98,bottom=.09,top=.92,hspace=.16,wspace=.30)
        subject_stats=[]
        for col,(name,(y,unit,note)) in enumerate(outputs.items()):
            Y=np.array(Ys[name]);H,Z,snr,snr_adj,coh,p,padj,psd=metrics(U,Y)
            for cyc in range(cycles):
                axes[0,col].loglog(f,abs(Z[cyc]),'o-',alpha=.25,lw=.6,color='gray')
                axes[1,col].semilogx(f,np.angle(Z[cyc],deg=True),'o',alpha=.35,color='gray',ms=3)
                for j,hz in enumerate(f):reps.append(dict(participant=alias,trial_id=tid,output=name,cycle=cyc+1,start_s=start+cyc*T,end_s=start+(cyc+1)*T,frequency_hz=hz,U_real=U[cyc,j].real,U_imag=U[cyc,j].imag,Y_real=Y[cyc,j].real,Y_imag=Y[cyc,j].imag,H_real=Z[cyc,j].real,H_imag=Z[cyc,j].imag,gain=abs(Z[cyc,j]),phase_deg=np.angle(Z[cyc,j],deg=True)))
            axes[0,col].loglog(f,abs(H),'o-',color='#0072B2',label='complex aggregate; gray: cycles');axes[0,col].set_title('COM acceleration Fy/M' if col==0 else 'Head IMU Z (raw counts)',fontsize=10)
            axes[1,col].semilogx(f,np.angle(H,deg=True),'o-',color='#0072B2');axes[1,col].set_ylim(-190,190)
            axes[2,col].semilogx(f,coh,'o-',color='#009E73');axes[2,col].set_ylim(0,1)
            axes[3,col].semilogx(f,snr,'o-',color='#D55E00',label='raw coherent / cycle variance');axes[3,col].semilogx(f,snr_adj,'x--',color='#8052a1',label='finite-K bias adjusted')
            for row in range(4):
                axes[row,col].grid(alpha=.25,which='both');axes[row,col].set_xticks(f,labels=[str(x) for x in f],rotation=65,fontsize=7)
            axes[0,col].set_ylabel(f'Gain [{unit}/(m/s)]');axes[1,col].set_ylabel('Phase (deg)');axes[2,col].set_ylabel(f'Coherence (K={cycles})');axes[3,col].set_ylabel('S/N (dB)');axes[3,col].set_xlabel('Frequency (Hz)')
            if col==0:axes[0,col].legend(fontsize=6);axes[3,col].legend(fontsize=6)
            for j,hz in enumerate(f):results.append(dict(participant=alias,trial_id=tid,output=name,frequency_hz=hz,gain=abs(H[j]),gain_unit=unit+'/(m/s)',phase_deg=np.angle(H[j],deg=True),H_real=H[j].real,H_imag=H[j].imag,cycles=cycles,independent_trials=1,coherence=coh[j],snr_db=snr[j],snr_bias_adjusted_db=snr_adj[j],phase_circular_sd_deg=psd[j],nominal_null_p=p[j],BH_q_per_output=padj[j],mean_output_amplitude=abs(Y[:,j].mean()),output_unit=unit,axis_and_model_status=note))
            subject_stats.append(dict(output=name,median_snr_db=float(np.median(snr)),median_coherence=float(np.median(coh)),q_below_0p05_hz=f[padj<.05].tolist(),phase_sd_deg=np.round(psd,1).tolist()))
        fig.suptitle(alias+' | 1 trial, '+str(cycles)+' consecutive cycles; no time-shift optimization; COM position unavailable',fontsize=11)
        fig.savefig(out/(alias.replace(' ','_')+'_frequency_response.png'),dpi=150);plt.close(fig)
        ax=ax_input[index,0];ax.plot(qe,q['ap_velocity_m_s'],lw=.65,color='gray');
        for cyc in range(cycles):ax.axvspan(start+cyc*T,start+(cyc+1)*T,color=['#0072B2','#009E73'][cyc%2],alpha=.15)
        ax.set(title=alias+' | shading: four complete periods',ylabel='Background AP velocity (m/s)',xlim=(0,c['duration_s']));ax.grid(alpha=.2)
        inv.update(status='included',analysis_start_s=start,analysis_end_s=stop,cycles=cycles,period_s=T,reason='complete; conditional sensor AP signs, unvalidated static-Fz mass')
        clock=m.get('quest_clock') or {}
        sync.append(dict(participant=alias,trial_id=tid,quest_frames=len(qe),force_samples=len(ft),quest_median_dt_s=float(np.median(np.diff(qt))),quest_max_gap_s=float(np.max(np.diff(qt))),force_median_dt_s=float(np.median(np.diff(ft))),force_max_gap_s=float(np.max(np.diff(ft))),clock_drift_ppm=clock.get('drift_ppm'),clock_selected_residual_rms_s=clock.get('selected_residual_rms_s'),mass_kg=mass,mass_method=m['weight_estimate']['method'],height_m=height,alignment='stored Quest PC mapping + estimated force sample times; no shift',physical_latency='uncalibrated'))
        for pth in [mp,qp,fp]:provenance.append(dict(participant=alias,file=str(pth.relative_to(root)),sha256=hashlib.sha256(pth.read_bytes()).hexdigest()))
        subject_summary.append(dict(participant=alias,frequency_hz=f.tolist(),period_s=T,cycles=cycles,outputs=subject_stats))
      except Exception as exc:
        inv.update(status='excluded',reason=str(exc));failures.append(dict(metadata=str(mp),error=str(exc)))
    ax_input[-1,0].set_xlabel('Stimulus elapsed time (s)');fig_input.savefig(out/'stimulus_and_intervals.png',dpi=150);plt.close(fig_input)
    for filename,rows in [('trial_inventory.csv',inventory),('signal_schema.csv',schema),('frequency_response.csv',results),('cycle_coefficients.csv',reps),('input_log_checks.csv',inputs),('synchronization_audit.csv',sync),('source_hashes.csv',provenance)]:write(out/filename,rows)
    (out/'summary.json').write_text(json.dumps(dict(subjects=subject_summary,failures=failures,canonical_main_trials=len(selected),analyzed_trials=sum(x['status']=='included' for x in inventory)),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(output=str(out),trials=len(selected),successful=len(subject_summary),failures=failures,rows=len(results)),indent=2))
    if failures:raise RuntimeError('Some target trials failed; see summary')

if __name__=='__main__':main()
