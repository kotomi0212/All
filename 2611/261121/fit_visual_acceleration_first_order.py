"""Identify a direct first-order visual-acceleration -> COM-acceleration model."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "visual_acceleration_first_order_fit"
FREQUENCIES = np.arange(0.05, 1.5000001, 0.05)
FIT_MAX_HZ = 0.70
FIT_MASK = FREQUENCIES <= FIT_MAX_HZ + 1e-12


def complex_amplitude(t: np.ndarray, y: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
    t = np.asarray(t, float) - float(t[0])
    y = np.asarray(y, float)
    trend = np.column_stack((np.ones_like(t), t))
    y = y - trend @ np.linalg.lstsq(trend, y, rcond=None)[0]
    result = []
    for frequency in frequencies:
        omega = 2 * np.pi * frequency
        basis = np.column_stack((np.cos(omega * t), np.sin(omega * t)))
        cosine, sine = np.linalg.lstsq(basis, y, rcond=None)[0]
        result.append(cosine - 1j * sine)
    return np.asarray(result)


def model(frequencies: np.ndarray, tv: float, gamma: float) -> np.ndarray:
    s = 1j * 2 * np.pi * frequencies
    return -gamma / (tv * s + 1)


def phase_difference(predicted: np.ndarray, measured: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (np.angle(predicted) - np.angle(measured))))


def fit(measured: np.ndarray):
    frequencies = FREQUENCIES[FIT_MASK]
    observed = measured[FIT_MASK]

    def residual(theta: np.ndarray) -> np.ndarray:
        tv, gamma = np.exp(theta)
        predicted = model(frequencies, tv, gamma)
        gain_error = np.log(np.maximum(np.abs(predicted), 1e-15)) - np.log(
            np.maximum(np.abs(observed), 1e-15)
        )
        phase_error = phase_difference(predicted, observed)
        return np.r_[gain_error, phase_error]

    result = least_squares(
        residual,
        x0=np.log([0.25, 0.05]),
        bounds=(np.log([0.001, 1e-8]), np.log([1000.0, 100.0])),
        max_nfev=30000,
    )
    tv, gamma = np.exp(result.x)
    predicted = model(frequencies, tv, gamma)
    gain_rmse_db = float(
        np.sqrt(np.mean((20 * np.log10(np.abs(predicted) / np.abs(observed))) ** 2))
    )
    phase_rmse_deg = float(
        np.sqrt(np.mean(np.degrees(phase_difference(predicted, observed)) ** 2))
    )
    return float(tv), float(gamma), gain_rmse_db, phase_rmse_deg


def load_trial(folder: Path):
    stem = folder.name
    quest = pd.read_csv(folder / f"{stem}_quest.csv")
    force = pd.read_csv(folder / f"{stem}_combalance.csv")
    total_mass = float(quest["body_mass_kg"].iloc[0])

    quest_time = quest["pc_mapped_time_s"].to_numpy(float)
    quest_elapsed = quest["stimulus_elapsed_s"].to_numpy(float)
    start_pc, end_pc = np.interp([10.0, 90.0], quest_elapsed, quest_time)
    force_time = force["pc_estimated_sample_time_s"].to_numpy(float)
    keep = (force_time >= start_pc) & (force_time < end_pc)
    time = force_time[keep] - start_pc

    visual_acceleration = np.interp(
        force_time[keep],
        quest_time,
        quest["ap_acceleration_m_s2"].to_numpy(float),
    )
    # 90-degree CCW participant orientation: AP force is plate Fy.
    com_acceleration = force.loc[keep, "sensor_Fy_N"].to_numpy(float) / total_mass
    input_amplitude = complex_amplitude(time, visual_acceleration, FREQUENCIES)
    output_amplitude = complex_amplitude(time, com_acceleration, FREQUENCIES)
    return {
        "trial": stem,
        "subject": str(quest["subject"].iloc[0]),
        "trial_number": int(quest["trial_number"].iloc[0]),
        "frf": output_amplitude / input_amplitude,
    }


def main():
    OUT.mkdir(exist_ok=True)
    trials = [
        load_trial(folder)
        for folder in sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
    ]
    dense = np.linspace(0.05, 1.50, 1000)
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True)
    summaries = []
    frequency_rows = []

    for column, trial in enumerate(trials):
        observed = trial["frf"]
        tv, gamma, gain_rmse_db, phase_rmse_deg = fit(observed)
        dense_model = model(dense, tv, gamma)
        point_model = model(FREQUENCIES, tv, gamma)
        gain_ax, phase_ax = axes[:, column]

        for ax in (gain_ax, phase_ax):
            ax.axvspan(0.05, FIT_MAX_HZ, color="#e8f3ec", alpha=0.8, zorder=0)
            ax.axvspan(FIT_MAX_HZ, 1.50, color="#f4f4f4", alpha=0.8, zorder=0)
            ax.axvline(FIT_MAX_HZ, color="#555555", linestyle="--", linewidth=1)
            ax.grid(True, which="both", alpha=0.25)
            ax.set_xlim(0.03, 1.52)

        gain_ax.plot(dense, np.abs(dense_model), color="#d62728", linewidth=2, label="model")
        gain_ax.scatter(FREQUENCIES[FIT_MASK], np.abs(observed[FIT_MASK]), color="#1f77b4", s=28, label="measured (fit)", zorder=3)
        gain_ax.scatter(FREQUENCIES[~FIT_MASK], np.abs(observed[~FIT_MASK]), facecolors="none", edgecolors="#777777", s=22, label="measured (>0.70 Hz)", zorder=3)
        gain_ax.set_yscale("log")
        gain_ax.set_title(f"{trial['subject']} trial {trial['trial_number']}")
        gain_ax.text(
            0.97, 0.96,
            f"$T_v$={tv:.3g} s\n$\\gamma_v$={gamma:.3g}\ngain RMSE={gain_rmse_db:.1f} dB\nphase RMSE={phase_rmse_deg:.1f}°",
            transform=gain_ax.transAxes, ha="right", va="top", fontsize=9,
            bbox={"facecolor":"white","alpha":0.82,"edgecolor":"none"},
        )

        phase_ax.plot(dense, np.degrees(np.angle(dense_model)), color="#d62728", linewidth=2)
        phase_ax.scatter(FREQUENCIES[FIT_MASK], np.degrees(np.angle(observed[FIT_MASK])), color="#1f77b4", s=28, zorder=3)
        phase_ax.scatter(FREQUENCIES[~FIT_MASK], np.degrees(np.angle(observed[~FIT_MASK])), facecolors="none", edgecolors="#777777", s=22, zorder=3)
        phase_ax.set_ylim(-190, 190)
        phase_ax.set_yticks([-180, -90, 0, 90, 180])
        phase_ax.set_xlabel("Frequency [Hz]")
        if column == 0:
            gain_ax.set_ylabel("Acceleration gain [-]")
            phase_ax.set_ylabel("Phase [deg]")

        summaries.append({
            "trial": trial["trial"], "subject": trial["subject"],
            "trial_number": trial["trial_number"], "T_v_s": tv,
            "gamma_v": gamma, "W_v": gamma / tv,
            "gain_rmse_db_0p05_0p70": gain_rmse_db,
            "phase_rmse_deg_0p05_0p70": phase_rmse_deg,
        })
        for frequency, measured_value, model_value in zip(FREQUENCIES, observed, point_model):
            frequency_rows.append({
                "trial": trial["trial"], "frequency_hz": frequency,
                "used_for_fit": frequency <= FIT_MAX_HZ + 1e-12,
                "measured_gain": abs(measured_value), "model_gain": abs(model_value),
                "measured_phase_deg": np.degrees(np.angle(measured_value)),
                "model_phase_deg": np.degrees(np.angle(model_value)),
            })

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.96), ncol=3, frameon=False)
    fig.suptitle("First-order model: visual acceleration to direct COM acceleration", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(OUT / "visual_to_com_acceleration_gain_phase.png", dpi=200)
    plt.close(fig)

    pd.DataFrame(summaries).to_csv(OUT / "identified_parameters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(frequency_rows).to_csv(OUT / "frequency_response_values.csv", index=False, encoding="utf-8-sig")
    (OUT / "model_equation.md").write_text(
        "# Model\n\n"
        "$$T_v\\dot a_{COM}(t)+a_{COM}(t)=-\\gamma_v a_{visual}(t)$$\n\n"
        "$$G(s)=\\frac{A_{COM}(s)}{A_{visual}(s)}=-\\frac{\\gamma_v}{T_vs+1}$$\n\n"
        "Measured output: $a_{COM}=F_y/M$. Identification uses gain and phase at 0.05-0.70 Hz.\n",
        encoding="utf-8",
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
