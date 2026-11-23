"""Compare plausible sagittal force-plate axis mappings for the 261121 trials."""

from pathlib import Path
import numpy as np
import pandas as pd

from analyze_visual_identification import (
    ROOT, FREQUENCIES_HZ, complex_amplitude, fit_visual_block, model_frf,
)


MAPPINGS = {
    "forward=+Y (Rx=+Fy, Ny=-Mx)": ("sensor_Fy_N", 1.0, "sensor_Mx_Nm", -1.0),
    "forward=-Y (Rx=-Fy, Ny=+Mx)": ("sensor_Fy_N", -1.0, "sensor_Mx_Nm", 1.0),
    "forward=+X (Rx=+Fx, Ny=+My)": ("sensor_Fx_N", 1.0, "sensor_My_Nm", 1.0),
    "forward=-X (Rx=-Fx, Ny=-My)": ("sensor_Fx_N", -1.0, "sensor_My_Nm", -1.0),
}


def trial_frf(folder: Path, mapping):
    stem = folder.name
    q = pd.read_csv(folder / f"{stem}_quest.csv")
    f = pd.read_csv(folder / f"{stem}_combalance.csv")
    H = float(q.height_cm.iloc[0]) / 100
    M = float(q.body_mass_kg.iloc[0])
    mb, jb, lb, Lf = .978*M, .0425*M*H**2, .531*H, .038*H
    A = jb/lb + mb*lb
    C = jb/(mb*lb) + lb + Lf
    qt, qe = q.pc_mapped_time_s.to_numpy(), q.stimulus_elapsed_s.to_numpy()
    start, end = np.interp([10., 90.], qe, qt)
    ft = f.pc_estimated_sample_time_s.to_numpy()
    keep = (ft >= start) & (ft < end)
    t = ft[keep] - start
    v = np.interp(ft[keep], qt, q.ap_velocity_m_s)
    fc, fs, mc, ms = mapping
    rx = fs * f.loc[keep, fc].to_numpy()
    ny = ms * f.loc[keep, mc].to_numpy()
    xb = (ny-C*rx)/mb
    frf = complex_amplitude(t, xb, FREQUENCIES_HZ) / complex_amplitude(t, v, FREQUENCIES_HZ)
    return str(q.subject.iloc[0]), frf, {"mb": mb, "lb": lb, "inertia_term": A}


rows=[]
folders=sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
for label,mapping in MAPPINGS.items():
    values={}
    pars={}
    for folder in folders:
        subject, frf, p=trial_frf(folder,mapping)
        values.setdefault(subject,[]).append(frf)
        pars[subject]=p
    for subject,reps in values.items():
        a,b=reps
        mean=(a+b)/2
        repeat_nrmse=np.sqrt(np.sum(np.abs(a-b)**2)/np.sum(np.abs(mean)**2))
        tv,gamma,fit=fit_visual_block(FREQUENCIES_HZ,mean,pars[subject])
        rows.append({"mapping":label,"subject":subject,"repeat_difference_nrmse":repeat_nrmse,
                     "T_v_s":tv,"gamma_v_s":gamma,"fit_percent":fit})

out=pd.DataFrame(rows)
out.to_csv(ROOT/'visual_identification_results'/'axis_diagnostics.csv',index=False,encoding='utf-8-sig')
print(out.to_string(index=False))
