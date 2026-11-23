"""Check both force-plate horizontal axes as direct COM acceleration outputs."""
from pathlib import Path
import numpy as np
import pandas as pd
from analyze_visual_identification import ROOT, FREQUENCIES_HZ, complex_amplitude, fit_visual_block

rows=[]
for axis in ("sensor_Fx_N", "sensor_Fy_N"):
    for sign in (1., -1.):
        for folder in sorted(ROOT.glob("20261121_*_90s_multicosine_*")):
            stem=folder.name
            q=pd.read_csv(folder/f"{stem}_quest.csv")
            f=pd.read_csv(folder/f"{stem}_combalance.csv")
            H=float(q.height_cm.iloc[0])/100; M=float(q.body_mass_kg.iloc[0])
            mb=.978*M; jb=.0425*M*H**2; lb=.531*H; A=jb/lb+mb*lb
            qt=q.pc_mapped_time_s.to_numpy(); qe=q.stimulus_elapsed_s.to_numpy()
            start,end=np.interp([10.,90.],qe,qt); ft=f.pc_estimated_sample_time_s.to_numpy()
            keep=(ft>=start)&(ft<end); t=ft[keep]-start
            v=np.interp(ft[keep],qt,q.ap_velocity_m_s)
            acc=sign*f.loc[keep,axis].to_numpy()/M
            frf=complex_amplitude(t,acc,FREQUENCIES_HZ)/complex_amplitude(t,v,FREQUENCIES_HZ)
            tv,gamma,fit=fit_visual_block(FREQUENCIES_HZ,frf,{"mb":mb,"lb":lb,"inertia_term":A})
            rows.append({"axis":axis,"sign":sign,"trial":stem,"T_v_s":tv,"gamma_v_s":gamma,"fit_percent":fit})
out=pd.DataFrame(rows)
path=ROOT/'visual_identification_results_acceleration'/'axis_diagnostics.csv'
out.to_csv(path,index=False,encoding='utf-8-sig')
print(out.to_string(index=False))
