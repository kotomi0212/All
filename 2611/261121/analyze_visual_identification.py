"""Identify T_v and gamma_v from the 2026-11-21 VR/force-plate trials.

The known balance controller follows Sonobe's note (2026-11-12):
Kp=10.53, Kd=3.97, td=0.155 s.  The visual block is
    T_v * d(xbar)/dt + xbar = -gamma_v * v_vr.

The participant faced 90 deg counter-clockwise from the previous force-plate
orientation.  Assuming the previous plate +X was participant-forward, the
participant sagittal components are R_x=plate F_y and N_y=-plate M_x.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
FREQUENCIES_HZ = np.arange(0.05, 1.5000001, 0.05)
KP_NORM = 10.53
KD_NORM = 3.97
DELAY_S = 0.155
G = 9.80665
OUTPUT_MODE = "velocity"
FIT_MODE = "gain_only"


def complex_amplitude(t: np.ndarray, y: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    """Least-squares complex sinusoid amplitude after removing mean/trend."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    t = t - t[0]
    # Remove offset and linear drift without distorting the exact excited bins.
    trend = np.column_stack((np.ones_like(t), t))
    y = y - trend @ np.linalg.lstsq(trend, y, rcond=None)[0]
    out = []
    for f in frequencies:
        w = 2 * np.pi * f
        basis = np.column_stack((np.cos(w * t), np.sin(w * t)))
        a, b = np.linalg.lstsq(basis, y, rcond=None)[0]
        out.append(a - 1j * b)
    return np.asarray(out)


def model_frf(freq_hz: np.ndarray, tv: float, gamma: float, *, mb: float, lb: float, inertia_term: float) -> np.ndarray:
    """FRF from VR velocity to the selected COM output."""
    s = 1j * 2 * np.pi * freq_hz
    kp = KP_NORM * inertia_term
    kd = KD_NORM * inertia_term
    delay = np.exp(-DELAY_S * s)
    # Requested reduced model: omit (J_b/l_b + m_b*l_b) * x_b_ddot.
    # The support-acceleration term is also absent because the platform is fixed.
    denominator = (kp + kd * s) * delay - mb * G
    displacement_frf = -kp * gamma * delay / ((tv * s + 1) * denominator)
    if OUTPUT_MODE == "acceleration":
        return s**2 * displacement_frf
    if OUTPUT_MODE == "velocity":
        return s * displacement_frf
    return displacement_frf


def fit_visual_block(freq: np.ndarray, measured: np.ndarray, params: dict[str, float]) -> tuple[float, float, float]:
    """Fit Tv and gamma to FRF gain only; phase is intentionally discarded."""
    # For acceleration output, use physical-unit residuals.  This intentionally
    # gives the larger high-frequency acceleration response its natural weight.
    if OUTPUT_MODE in {"velocity", "acceleration"}:
        scale = np.ones_like(np.abs(measured))
    else:
        scale = np.maximum(np.abs(measured), np.nanmedian(np.abs(measured)) * 0.25)

    def residual(theta: np.ndarray) -> np.ndarray:
        tv = math.exp(theta[0])
        gamma = math.exp(theta[1])
        predicted = model_frf(freq, tv, gamma, **params)
        return (np.abs(predicted) - np.abs(measured)) / scale

    result = least_squares(
        residual,
        x0=np.array([math.log(2.5), math.log(0.04)]),
        bounds=([math.log(0.001), math.log(1e-8)], [math.log(1000.0), math.log(100.0)]),
        max_nfev=20000,
    )
    tv = math.exp(result.x[0])
    gamma = math.exp(result.x[1])
    predicted = model_frf(freq, tv, gamma, **params)
    measured_gain = np.abs(measured)
    predicted_gain = np.abs(predicted)
    ss_error = float(np.sum((measured_gain - predicted_gain) ** 2))
    ss_total = float(np.sum((measured_gain - np.mean(measured_gain)) ** 2))
    fit_percent = 100 * (1 - ss_error / ss_total) if ss_total > 0 else float("nan")
    return tv, gamma, fit_percent


def analyze_trial(folder: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    stem = folder.name
    quest = pd.read_csv(folder / f"{stem}_quest.csv")
    force = pd.read_csv(folder / f"{stem}_combalance.csv")
    metadata = json.loads((folder / f"{stem}_metadata.json").read_text(encoding="utf-8"))

    height_m = float(quest["height_cm"].iloc[0]) / 100
    total_mass = float(quest["body_mass_kg"].iloc[0])
    mb = 0.978 * total_mass
    jb = 0.0425 * total_mass * height_m**2
    lb = 0.531 * height_m
    lf = 0.015 * height_m
    Lf = 0.038 * height_m
    inertia_term = jb / lb + mb * lb
    com_force_coefficient = jb / (mb * lb) + lb + Lf

    # Align by the shared PC monotonic clock and retain the pure multicosine
    # stationary interval.  A tiny endpoint margin avoids partial samples.
    q_time = quest["pc_mapped_time_s"].to_numpy(float)
    q_elapsed = quest["stimulus_elapsed_s"].to_numpy(float)
    start_pc = float(np.interp(10.0, q_elapsed, q_time))
    end_pc = float(np.interp(90.0, q_elapsed, q_time))
    f_time = force["pc_estimated_sample_time_s"].to_numpy(float)
    keep = (f_time >= start_pc) & (f_time < end_pc)
    t_pc = f_time[keep]
    t = t_pc - start_pc

    v_vr = np.interp(t_pc, q_time, quest["ap_velocity_m_s"].to_numpy(float))
    rx = force.loc[keep, "sensor_Fy_N"].to_numpy(float)
    ny = -force.loc[keep, "sensor_Mx_Nm"].to_numpy(float)
    xb = (ny - com_force_coefficient * rx) / mb

    input_amp = complex_amplitude(t, v_vr, FREQUENCIES_HZ)
    if OUTPUT_MODE == "acceleration":
        # Newton's law gives COM acceleration directly from the horizontal
        # ground-reaction force.  This avoids differentiating the very small
        # force/moment-derived COM displacement twice.
        com_acceleration = rx / total_mass
        output_amp = complex_amplitude(t, com_acceleration, FREQUENCIES_HZ)
    elif OUTPUT_MODE == "velocity":
        # Integrate the directly observed COM acceleration in the frequency
        # domain. There is no DC division because all fitted bins are > 0 Hz.
        s = 1j * 2 * np.pi * FREQUENCIES_HZ
        com_acceleration = rx / total_mass
        acceleration_amp = complex_amplitude(t, com_acceleration, FREQUENCIES_HZ)
        output_amp = acceleration_amp / s
    else:
        output_amp = complex_amplitude(t, xb, FREQUENCIES_HZ)
    measured = output_amp / input_amp
    model_params = {"mb": mb, "lb": lb, "inertia_term": inertia_term}
    tv, gamma, fit_percent = fit_visual_block(FREQUENCIES_HZ, measured, model_params)
    predicted = model_frf(FREQUENCIES_HZ, tv, gamma, **model_params)

    summary = {
        "trial": stem,
        "subject": str(quest["subject"].iloc[0]),
        "trial_number": int(quest["trial_number"].iloc[0]),
        "height_m": height_m,
        "body_mass_kg": total_mass,
        "samples": int(keep.sum()),
        "duration_s": float(t[-1] - t[0]),
        "T_v_s": tv,
        "gamma_v_s": gamma,
        "W_v": gamma / tv,
        "fit_percent_gain": fit_percent,
        "axis_mapping": "Rx_subject=Fy_plate; Ny_subject=-Mx_plate",
        "fit_band_hz": "0.05-1.50",
        "output_mode": OUTPUT_MODE,
        "fit_mode": FIT_MODE,
    }
    rows = []
    for f, obs, pred, u in zip(FREQUENCIES_HZ, measured, predicted, input_amp):
        rows.append(
            {
                "trial": stem,
                "frequency_hz": f,
                "input_velocity_amplitude_m_s": abs(u),
                "measured_gain": abs(obs),
                "measured_phase_deg": np.degrees(np.angle(obs)),
                "model_gain": abs(pred),
                "model_phase_deg": np.degrees(np.angle(pred)),
                "complex_error_abs": abs(obs - pred),
            }
        )
    return summary, rows


def main() -> None:
    folders = sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
    summaries: list[dict[str, object]] = []
    frequency_rows: list[dict[str, object]] = []
    for folder in folders:
        summary, rows = analyze_trial(folder)
        summaries.append(summary)
        frequency_rows.extend(rows)

    output = ROOT / f"visual_identification_results_{OUTPUT_MODE}"
    output.mkdir(exist_ok=True)
    pd.DataFrame(summaries).to_csv(output / "trial_parameters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(frequency_rows).to_csv(output / "frequency_fits.csv", index=False, encoding="utf-8-sig")

    subject_summary = (
        pd.DataFrame(summaries)
        .groupby("subject")[["T_v_s", "gamma_v_s", "W_v", "fit_percent_gain"]]
        .agg(["mean", "std"])
    )
    subject_summary.to_csv(output / "subject_parameters.csv", encoding="utf-8-sig")
    print(pd.DataFrame(summaries).to_string(index=False))
    print("\nSubject summary:\n", subject_summary.to_string())


if __name__ == "__main__":
    main()
