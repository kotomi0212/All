"""Compare gain-only model fits for COM displacement, velocity, acceleration.

Only 0.05-0.70 Hz is used for identification. Measured points up to 1.50 Hz
are displayed to show the response to every component in the current stimulus.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "com_output_model_comparison"
ALL_FREQUENCIES = np.arange(0.05, 1.5000001, 0.05)
FIT_MAX_HZ = 0.70
FIT_MASK = ALL_FREQUENCIES <= FIT_MAX_HZ + 1e-12
KP_NORM = 10.53
KD_NORM = 3.97
DELAY_S = 0.155
GRAVITY = 9.80665
OUTPUTS = {
    "displacement": {"power": 0, "gain_unit": "s", "label": "COM displacement"},
    "velocity": {"power": 1, "gain_unit": "-", "label": "COM velocity"},
    "acceleration": {"power": 2, "gain_unit": "1/s", "label": "COM acceleration"},
}


def complex_amplitude(t: np.ndarray, y: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float) - float(t[0])
    y = np.asarray(y, float)
    trend = np.column_stack((np.ones_like(t), t))
    y = y - trend @ np.linalg.lstsq(trend, y, rcond=None)[0]
    amplitudes = []
    for frequency in frequencies:
        omega = 2 * np.pi * frequency
        basis = np.column_stack((np.cos(omega * t), np.sin(omega * t)))
        cosine, sine = np.linalg.lstsq(basis, y, rcond=None)[0]
        amplitudes.append(cosine - 1j * sine)
    return np.asarray(amplitudes)


def displacement_model(
    frequencies: np.ndarray,
    tv: float,
    gamma: float,
    *,
    body_mass: float,
    inertia_term: float,
) -> np.ndarray:
    """VR velocity -> COM displacement, with Eq. (1) inertia term omitted."""
    s = 1j * 2 * np.pi * frequencies
    kp = KP_NORM * inertia_term
    kd = KD_NORM * inertia_term
    delay = np.exp(-DELAY_S * s)
    denominator = (kp + kd * s) * delay - body_mass * GRAVITY
    return -kp * gamma * delay / ((tv * s + 1) * denominator)


def output_model(
    frequencies: np.ndarray,
    tv: float,
    gamma: float,
    output_power: int,
    parameters: dict[str, float],
) -> np.ndarray:
    s = 1j * 2 * np.pi * frequencies
    return s**output_power * displacement_model(frequencies, tv, gamma, **parameters)


def fit_gain(
    measured_frf: np.ndarray,
    output_power: int,
    parameters: dict[str, float],
) -> tuple[float, float, float]:
    measured_gain = np.abs(measured_frf[FIT_MASK])
    frequencies = ALL_FREQUENCIES[FIT_MASK]

    def residual(theta: np.ndarray) -> np.ndarray:
        tv, gamma = np.exp(theta)
        predicted = output_model(frequencies, tv, gamma, output_power, parameters)
        # Physical-unit, gain-only residual. Changing output therefore changes
        # frequency weighting, which is the comparison requested here.
        return np.abs(predicted) - measured_gain

    result = least_squares(
        residual,
        x0=np.log([0.25, 0.03]),
        bounds=(np.log([0.001, 1e-8]), np.log([1000.0, 100.0])),
        max_nfev=20000,
    )
    tv, gamma = np.exp(result.x)
    predicted_gain = np.abs(output_model(frequencies, tv, gamma, output_power, parameters))
    ss_error = float(np.sum((measured_gain - predicted_gain) ** 2))
    ss_total = float(np.sum((measured_gain - np.mean(measured_gain)) ** 2))
    fit_percent = 100 * (1 - ss_error / ss_total) if ss_total > 0 else float("nan")
    return float(tv), float(gamma), fit_percent


def load_trial(folder: Path) -> dict[str, object]:
    stem = folder.name
    quest = pd.read_csv(folder / f"{stem}_quest.csv")
    force = pd.read_csv(folder / f"{stem}_combalance.csv")
    height = float(quest["height_cm"].iloc[0]) / 100
    total_mass = float(quest["body_mass_kg"].iloc[0])
    body_mass = 0.978 * total_mass
    body_inertia = 0.0425 * total_mass * height**2
    com_height = 0.531 * height
    ankle_height = 0.038 * height
    inertia_term = body_inertia / com_height + body_mass * com_height
    force_coefficient = body_inertia / (body_mass * com_height) + com_height + ankle_height

    quest_time = quest["pc_mapped_time_s"].to_numpy(float)
    quest_elapsed = quest["stimulus_elapsed_s"].to_numpy(float)
    start_pc, end_pc = np.interp([10.0, 90.0], quest_elapsed, quest_time)
    force_time = force["pc_estimated_sample_time_s"].to_numpy(float)
    keep = (force_time >= start_pc) & (force_time < end_pc)
    time = force_time[keep] - start_pc
    vr_velocity = np.interp(
        force_time[keep], quest_time, quest["ap_velocity_m_s"].to_numpy(float)
    )

    # Participant turned 90 degrees CCW: sagittal Rx=plate Fy, Ny=-plate Mx.
    horizontal_force = force.loc[keep, "sensor_Fy_N"].to_numpy(float)
    sagittal_moment = -force.loc[keep, "sensor_Mx_Nm"].to_numpy(float)
    com_displacement = (
        sagittal_moment - force_coefficient * horizontal_force
    ) / body_mass

    input_amplitude = complex_amplitude(time, vr_velocity, ALL_FREQUENCIES)
    displacement_amplitude = complex_amplitude(time, com_displacement, ALL_FREQUENCIES)
    displacement_frf = displacement_amplitude / input_amplitude
    return {
        "trial": stem,
        "subject": str(quest["subject"].iloc[0]),
        "trial_number": int(quest["trial_number"].iloc[0]),
        "displacement_frf": displacement_frf,
        "parameters": {"body_mass": body_mass, "inertia_term": inertia_term},
    }


def write_method_note() -> None:
    note = r"""# COM出力別モデル比較

## 同定条件

- 入力: VRのAP速度 $V_{VR}(s)$
- 同定周波数: 0.05-0.70 Hz（14点）
- 0.75-1.50 Hz: グラフには表示するが同定には不使用
- 位相: 使用しない（ゲインのみ）
- 既知値: $K_p=10.53$, $K_d=3.97$, $t_d=0.155$ s
- 式(1)の身体慣性項 $(J_b/l_b+m_bl_b)\ddot{x}_b$ は指示どおり削除

## 基本となる変位出力

式(2), (3)から、VR速度からCOM変位への伝達関数を

$$
G_x(s)=\frac{X_b(s)}{V_{VR}(s)}
=\frac{-k_p\gamma_v e^{-t_ds}}
{(T_vs+1)\{(k_p+k_ds)e^{-t_ds}-m_bg\}}
$$

とした。

## 速度出力への変更

出力を $X_b(s)$ から $\dot{X}_b(s)=sX_b(s)$ に変更した。したがって

$$
G_v(s)=\frac{\dot{X}_b(s)}{V_{VR}(s)}=sG_x(s)
$$

実測値も、COM変位の各周波数の複素フーリエ係数に $s=j2\pi f$ を掛けて速度へ変換した。

## 加速度出力への変更

出力を $\ddot{X}_b(s)=s^2X_b(s)$ に変更した。したがって

$$
G_a(s)=\frac{\ddot{X}_b(s)}{V_{VR}(s)}=s^2G_x(s)
$$

実測値も、COM変位の各周波数の複素フーリエ係数に $(j2\pi f)^2$ を掛けて加速度へ変換した。

各出力について、0.05-0.70 Hzの実測ゲインとモデルゲインの物理単位での二乗誤差を最小化し、$T_v,\gamma_v$を別々に同定した。
"""
    (OUT / "method_and_equations.md").write_text(note, encoding="utf-8")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    folders = sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
    trials = [load_trial(folder) for folder in folders]
    result_rows: list[dict[str, object]] = []
    frequency_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=True)
    dense_frequencies = np.linspace(0.05, 1.50, 800)

    for row_index, (output_name, specification) in enumerate(OUTPUTS.items()):
        power = int(specification["power"])
        for column_index, trial in enumerate(trials):
            ax = axes[row_index, column_index]
            s = 1j * 2 * np.pi * ALL_FREQUENCIES
            measured_frf = s**power * trial["displacement_frf"]
            tv, gamma, fit_percent = fit_gain(
                measured_frf, power, trial["parameters"]
            )
            predicted_dense = output_model(
                dense_frequencies, tv, gamma, power, trial["parameters"]
            )
            predicted_points = output_model(
                ALL_FREQUENCIES, tv, gamma, power, trial["parameters"]
            )

            ax.axvspan(0.05, FIT_MAX_HZ, color="#e8f3ec", alpha=0.75, zorder=0)
            ax.axvspan(FIT_MAX_HZ, 1.50, color="#f4f4f4", alpha=0.75, zorder=0)
            ax.plot(
                dense_frequencies,
                np.abs(predicted_dense),
                color="#d62728",
                linewidth=2,
                label="model",
            )
            ax.scatter(
                ALL_FREQUENCIES[FIT_MASK],
                np.abs(measured_frf[FIT_MASK]),
                s=28,
                color="#1f77b4",
                label="measured (fit)",
                zorder=3,
            )
            ax.scatter(
                ALL_FREQUENCIES[~FIT_MASK],
                np.abs(measured_frf[~FIT_MASK]),
                s=22,
                facecolors="none",
                edgecolors="#777777",
                label="measured (>0.70 Hz, not fitted)",
                zorder=3,
            )
            ax.axvline(FIT_MAX_HZ, color="#555555", linestyle="--", linewidth=1)
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25)
            ax.set_xlim(0.03, 1.52)
            ax.text(
                0.97,
                0.96,
                f"$T_v$={tv:.3g} s\n$\\gamma_v$={gamma:.3g}\nfit={fit_percent:.1f}%",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            )
            if row_index == 0:
                ax.set_title(f"{trial['subject']} trial {trial['trial_number']}")
            if column_index == 0:
                ax.set_ylabel(
                    f"{specification['label']} gain [{specification['gain_unit']}]"
                )
            if row_index == 2:
                ax.set_xlabel("Frequency [Hz]")

            result_rows.append(
                {
                    "output": output_name,
                    "trial": trial["trial"],
                    "subject": trial["subject"],
                    "trial_number": trial["trial_number"],
                    "T_v_s": tv,
                    "gamma_v": gamma,
                    "W_v": gamma / tv,
                    "fit_percent_gain_0p05_to_0p70_Hz": fit_percent,
                }
            )
            for frequency, observed, predicted in zip(
                ALL_FREQUENCIES, measured_frf, predicted_points
            ):
                frequency_rows.append(
                    {
                        "output": output_name,
                        "trial": trial["trial"],
                        "frequency_hz": frequency,
                        "used_for_fit": frequency <= FIT_MAX_HZ + 1e-12,
                        "measured_gain": abs(observed),
                        "model_gain": abs(predicted),
                    }
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "Gain-only model fits for three COM outputs (identification band: 0.05-0.70 Hz)",
        y=0.995,
        fontsize=15,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    fig.savefig(OUT / "model_vs_measured_three_outputs.png", dpi=200)
    plt.close(fig)

    results = pd.DataFrame(result_rows)
    results.to_csv(OUT / "identified_parameters_three_outputs.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(frequency_rows).to_csv(
        OUT / "model_vs_measured_values.csv", index=False, encoding="utf-8-sig"
    )
    write_method_note()
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
