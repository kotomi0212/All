"""Fit the visual model to directly measured COM acceleration, including phase."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "direct_com_acceleration_fit"
FREQUENCIES = np.arange(0.05, 1.5000001, 0.05)
FIT_MAX_HZ = 0.70
FIT_MASK = FREQUENCIES <= FIT_MAX_HZ + 1e-12
KP_NORM = 10.53
KD_NORM = 3.97
DELAY_S = 0.155
G = 9.80665


def complex_amplitude(t: np.ndarray, y: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float) - float(t[0])
    y = np.asarray(y, float)
    trend = np.column_stack((np.ones_like(t), t))
    y = y - trend @ np.linalg.lstsq(trend, y, rcond=None)[0]
    values = []
    for frequency in frequencies:
        omega = 2 * np.pi * frequency
        basis = np.column_stack((np.cos(omega * t), np.sin(omega * t)))
        cosine, sine = np.linalg.lstsq(basis, y, rcond=None)[0]
        values.append(cosine - 1j * sine)
    return np.asarray(values)


def acceleration_model(
    frequencies: np.ndarray,
    tv: float,
    gamma: float,
    *,
    body_mass: float,
    inertia_term: float,
) -> np.ndarray:
    """VR velocity -> COM acceleration; measured acceleration is not derived."""
    s = 1j * 2 * np.pi * frequencies
    kp = KP_NORM * inertia_term
    kd = KD_NORM * inertia_term
    delay = np.exp(-DELAY_S * s)
    # Requested reduced equation: omit body inertia A*x_b_ddot and platform
    # acceleration force. First obtain displacement response, then change only
    # the MODEL output to acceleration by multiplying its transfer function s^2.
    denominator = (kp + kd * s) * delay - body_mass * G
    displacement_response = (
        -kp * gamma * delay / ((tv * s + 1) * denominator)
    )
    return s**2 * displacement_response


def wrapped_phase_error(predicted: np.ndarray, measured: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (np.angle(predicted) - np.angle(measured))))


def fit_complex_frf(measured: np.ndarray, parameters: dict[str, float]):
    frequencies = FREQUENCIES[FIT_MASK]
    observed = measured[FIT_MASK]

    def residual(theta: np.ndarray) -> np.ndarray:
        tv, gamma = np.exp(theta)
        predicted = acceleration_model(frequencies, tv, gamma, **parameters)
        # Gain error is logarithmic; phase error is circular and in radians.
        # Both are dimensionless, preventing large high-frequency gains from
        # overwhelming the phase information.
        gain_error = np.log(np.maximum(np.abs(predicted), 1e-15)) - np.log(
            np.maximum(np.abs(observed), 1e-15)
        )
        phase_error = wrapped_phase_error(predicted, observed)
        return np.r_[gain_error, phase_error]

    result = least_squares(
        residual,
        x0=np.log([0.25, 0.03]),
        bounds=(np.log([0.001, 1e-8]), np.log([1000.0, 100.0])),
        max_nfev=30000,
    )
    tv, gamma = np.exp(result.x)
    predicted = acceleration_model(frequencies, tv, gamma, **parameters)
    gain_rmse_db = float(
        np.sqrt(np.mean((20 * np.log10(np.abs(predicted) / np.abs(observed))) ** 2))
    )
    phase_error_deg = np.degrees(wrapped_phase_error(predicted, observed))
    phase_rmse_deg = float(np.sqrt(np.mean(phase_error_deg**2)))
    return float(tv), float(gamma), gain_rmse_db, phase_rmse_deg


def load_trial(folder: Path):
    stem = folder.name
    quest = pd.read_csv(folder / f"{stem}_quest.csv")
    force = pd.read_csv(folder / f"{stem}_combalance.csv")
    height = float(quest["height_cm"].iloc[0]) / 100
    total_mass = float(quest["body_mass_kg"].iloc[0])
    body_mass = 0.978 * total_mass
    body_inertia = 0.0425 * total_mass * height**2
    com_height = 0.531 * height
    inertia_term = body_inertia / com_height + body_mass * com_height

    quest_time = quest["pc_mapped_time_s"].to_numpy(float)
    quest_elapsed = quest["stimulus_elapsed_s"].to_numpy(float)
    start_pc, end_pc = np.interp([10.0, 90.0], quest_elapsed, quest_time)
    force_time = force["pc_estimated_sample_time_s"].to_numpy(float)
    keep = (force_time >= start_pc) & (force_time < end_pc)
    time = force_time[keep] - start_pc
    vr_velocity = np.interp(
        force_time[keep], quest_time, quest["ap_velocity_m_s"].to_numpy(float)
    )
    # After the 90-degree CCW turn, participant AP force is plate Fy.
    com_acceleration = force.loc[keep, "sensor_Fy_N"].to_numpy(float) / total_mass
    input_amplitude = complex_amplitude(time, vr_velocity, FREQUENCIES)
    output_amplitude = complex_amplitude(time, com_acceleration, FREQUENCIES)
    return {
        "trial": stem,
        "subject": str(quest["subject"].iloc[0]),
        "trial_number": int(quest["trial_number"].iloc[0]),
        "measured_frf": output_amplitude / input_amplitude,
        "parameters": {"body_mass": body_mass, "inertia_term": inertia_term},
    }


def write_equation_note():
    text = r"""# 直接計測したCOM加速度による同定

## 実測出力

COM変位を微分せず、被験者前後方向の水平床反力から直接

$$
a_{COM}(t)=\frac{F_y(t)}{M}
$$

とした。

## モデル出力の変更

VR速度からCOM変位へのモデルを $G_x(s)=X_b(s)/V_{VR}(s)$ とすると、
モデルの出力だけを加速度へ変更して

$$
G_a(s)=\frac{A_b(s)}{V_{VR}(s)}=s^2G_x(s)
$$

とした。実測加速度に $s^2$ を掛ける操作は行っていない。

具体的には

$$
G_a(s)=
\frac{-s^2 k_p\gamma_v e^{-t_ds}}
{(T_vs+1)\{(k_p+k_ds)e^{-t_ds}-m_bg\}}
$$

である。0.05-0.70 Hzのゲインと位相を用いて $T_v,\gamma_v$ を同定した。

## 式(3)を加速度で読む場合

元の式

$$
T_v\dot{\bar{x}}_b+\bar{x}_b=-\gamma_v\dot{X}_s
$$

は、支持面速度からCOM平衡位置への一次遅れである。両辺を2回微分すると

$$
T_v\dot{\bar{a}}_b+\bar{a}_b=-\gamma_v\dddot{X}_s,
\qquad \bar{a}_b=\ddot{\bar{x}}_b
$$

となる。したがって、数学的に同じモデルを加速度出力で表すと、入力は支持面速度の
2階微分（支持面jerk）になる。単純な「視覚加速度→COM加速度」の一次遅れに置き換える
ことは、元式の変形ではなく、新しいモデル仮定になる。
"""
    (OUT / "acceleration_model_interpretation.md").write_text(text, encoding="utf-8")


def main():
    OUT.mkdir(exist_ok=True)
    folders = sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
    trials = [load_trial(folder) for folder in folders]
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    dense = np.linspace(0.05, 1.50, 1000)
    summary = []
    values = []

    for column, trial in enumerate(trials):
        measured = trial["measured_frf"]
        tv, gamma, gain_rmse_db, phase_rmse_deg = fit_complex_frf(
            measured, trial["parameters"]
        )
        model_dense = acceleration_model(dense, tv, gamma, **trial["parameters"])
        model_points = acceleration_model(FREQUENCIES, tv, gamma, **trial["parameters"])
        gain_ax, phase_ax = axes[:, column]
        for ax in (gain_ax, phase_ax):
            ax.axvspan(0.05, FIT_MAX_HZ, color="#e8f3ec", alpha=0.8, zorder=0)
            ax.axvspan(FIT_MAX_HZ, 1.50, color="#f4f4f4", alpha=0.8, zorder=0)
            ax.axvline(FIT_MAX_HZ, color="#555555", linestyle="--", linewidth=1)
            ax.grid(True, which="both", alpha=0.25)
            ax.set_xlim(0.03, 1.52)

        gain_ax.plot(dense, np.abs(model_dense), color="#d62728", linewidth=2, label="model")
        gain_ax.scatter(FREQUENCIES[FIT_MASK], np.abs(measured[FIT_MASK]), color="#1f77b4", s=28, label="measured (fit)", zorder=3)
        gain_ax.scatter(FREQUENCIES[~FIT_MASK], np.abs(measured[~FIT_MASK]), facecolors="none", edgecolors="#777777", s=22, label="measured (>0.70 Hz)", zorder=3)
        gain_ax.set_yscale("log")
        gain_ax.set_title(f"{trial['subject']} trial {trial['trial_number']}")
        gain_ax.text(0.97, 0.96, f"$T_v$={tv:.3g} s\n$\\gamma_v$={gamma:.3g}\ngain RMSE={gain_rmse_db:.1f} dB\nphase RMSE={phase_rmse_deg:.1f}°", transform=gain_ax.transAxes, ha="right", va="top", fontsize=9, bbox={"facecolor":"white","alpha":0.82,"edgecolor":"none"})

        measured_phase = np.degrees(np.angle(measured))
        model_phase_dense = np.degrees(np.angle(model_dense))
        phase_ax.plot(dense, model_phase_dense, color="#d62728", linewidth=2)
        phase_ax.scatter(FREQUENCIES[FIT_MASK], measured_phase[FIT_MASK], color="#1f77b4", s=28, zorder=3)
        phase_ax.scatter(FREQUENCIES[~FIT_MASK], measured_phase[~FIT_MASK], facecolors="none", edgecolors="#777777", s=22, zorder=3)
        phase_ax.set_ylim(-190, 190)
        phase_ax.set_yticks([-180, -90, 0, 90, 180])
        phase_ax.set_xlabel("Frequency [Hz]")

        if column == 0:
            gain_ax.set_ylabel("Acceleration gain [1/s]")
            phase_ax.set_ylabel("Phase [deg]")
        summary.append({"trial":trial["trial"],"subject":trial["subject"],"trial_number":trial["trial_number"],"T_v_s":tv,"gamma_v":gamma,"W_v":gamma/tv,"gain_rmse_db_0p05_0p70":gain_rmse_db,"phase_rmse_deg_0p05_0p70":phase_rmse_deg})
        for frequency, obs, pred in zip(FREQUENCIES, measured, model_points):
            values.append({"trial":trial["trial"],"frequency_hz":frequency,"used_for_fit":frequency <= FIT_MAX_HZ+1e-12,"measured_gain":abs(obs),"model_gain":abs(pred),"measured_phase_deg_wrapped":np.degrees(np.angle(obs)),"model_phase_deg_wrapped":np.degrees(np.angle(pred))})

    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5,0.96), ncol=3, frameon=False)
    fig.suptitle("Direct COM acceleration FRF: model vs measured (gain and phase; fit band 0.05-0.70 Hz)", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0,0,1,0.91))
    fig.savefig(OUT / "direct_com_acceleration_gain_phase.png", dpi=200)
    plt.close(fig)
    pd.DataFrame(summary).to_csv(OUT / "identified_parameters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(values).to_csv(OUT / "frequency_response_values.csv", index=False, encoding="utf-8-sig")
    write_equation_note()
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
