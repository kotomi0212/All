"""Gain-only fit of visual acceleration -> direct COM acceleration."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from fit_visual_acceleration_first_order import (
    ROOT,
    FREQUENCIES,
    FIT_MASK,
    FIT_MAX_HZ,
    load_trial,
    model,
)


OUT = ROOT / "visual_acceleration_first_order_gain_only"


def fit_gain_only(measured: np.ndarray):
    frequencies = FREQUENCIES[FIT_MASK]
    measured_gain = np.abs(measured[FIT_MASK])

    def residual(theta: np.ndarray) -> np.ndarray:
        tv, gamma = np.exp(theta)
        predicted_gain = np.abs(model(frequencies, tv, gamma))
        # Ordinary gain least squares; phase is not used. This is also the
        # residual definition used for the reported linear-gain R^2.
        return predicted_gain - measured_gain

    result = least_squares(
        residual,
        x0=np.log([0.25, 0.05]),
        bounds=(np.log([0.001, 1e-8]), np.log([1000.0, 100.0])),
        max_nfev=30000,
    )
    tv, gamma = np.exp(result.x)
    predicted_gain = np.abs(model(frequencies, tv, gamma))
    error = measured_gain - predicted_gain
    ss_error = float(np.sum(error**2))
    ss_total = float(np.sum((measured_gain - np.mean(measured_gain)) ** 2))
    r2_percent = 100 * (1 - ss_error / ss_total) if ss_total > 0 else float("nan")
    rmse_db = float(
        np.sqrt(np.mean((20 * np.log10(predicted_gain / measured_gain)) ** 2))
    )
    nrmse_percent = float(
        100 * np.sqrt(np.mean(error**2)) / np.mean(measured_gain)
    )
    return float(tv), float(gamma), r2_percent, rmse_db, nrmse_percent


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
        tv, gamma, r2_percent, rmse_db, nrmse_percent = fit_gain_only(observed)
        dense_model = model(dense, tv, gamma)
        point_model = model(FREQUENCIES, tv, gamma)
        gain_ax, phase_ax = axes[:, column]

        for ax in (gain_ax, phase_ax):
            ax.axvspan(0.05, FIT_MAX_HZ, color="#e8f3ec", alpha=0.8, zorder=0)
            ax.axvspan(FIT_MAX_HZ, 1.50, color="#f4f4f4", alpha=0.8, zorder=0)
            ax.axvline(FIT_MAX_HZ, color="#555555", linestyle="--", linewidth=1)
            ax.grid(True, which="both", alpha=0.25)
            ax.set_xlim(0.03, 1.52)

        gain_ax.plot(dense, np.abs(dense_model), color="#d62728", linewidth=2, label="model (gain-only fit)")
        gain_ax.scatter(FREQUENCIES[FIT_MASK], np.abs(observed[FIT_MASK]), color="#1f77b4", s=28, label="measured (fit)", zorder=3)
        gain_ax.scatter(FREQUENCIES[~FIT_MASK], np.abs(observed[~FIT_MASK]), facecolors="none", edgecolors="#777777", s=22, label="measured (>0.70 Hz)", zorder=3)
        gain_ax.set_yscale("log")
        gain_ax.set_title(f"{trial['subject']} trial {trial['trial_number']}")
        gain_ax.text(
            0.97,
            0.96,
            f"$T_v$={tv:.3g} s\n$\\gamma_v$={gamma:.3g}\n$R^2$={r2_percent:.1f}%\nRMSE={rmse_db:.1f} dB",
            transform=gain_ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )

        # Phase is displayed only as a diagnostic; it did not affect the fit.
        phase_ax.plot(dense, np.degrees(np.angle(dense_model)), color="#d62728", linewidth=2)
        phase_ax.scatter(FREQUENCIES[FIT_MASK], np.degrees(np.angle(observed[FIT_MASK])), color="#1f77b4", s=28, zorder=3)
        phase_ax.scatter(FREQUENCIES[~FIT_MASK], np.degrees(np.angle(observed[~FIT_MASK])), facecolors="none", edgecolors="#777777", s=22, zorder=3)
        phase_ax.set_ylim(-190, 190)
        phase_ax.set_yticks([-180, -90, 0, 90, 180])
        phase_ax.set_xlabel("Frequency [Hz]")
        if column == 0:
            gain_ax.set_ylabel("Acceleration gain [-]")
            phase_ax.set_ylabel("Phase [deg] (not fitted)")

        summaries.append(
            {
                "trial": trial["trial"],
                "subject": trial["subject"],
                "trial_number": trial["trial_number"],
                "T_v_s": tv,
                "gamma_v": gamma,
                "W_v": gamma / tv,
                "gain_R2_percent_0p05_0p70": r2_percent,
                "gain_RMSE_dB_0p05_0p70": rmse_db,
                "gain_NRMSE_percent_of_mean_0p05_0p70": nrmse_percent,
            }
        )
        for frequency, measured_value, model_value in zip(
            FREQUENCIES, observed, point_model
        ):
            frequency_rows.append(
                {
                    "trial": trial["trial"],
                    "frequency_hz": frequency,
                    "used_for_fit": frequency <= FIT_MAX_HZ + 1e-12,
                    "measured_gain": abs(measured_value),
                    "model_gain": abs(model_value),
                    "measured_phase_deg": np.degrees(np.angle(measured_value)),
                    "model_phase_deg_not_fitted": np.degrees(np.angle(model_value)),
                }
            )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.96),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        "Gain-only first-order fit: visual acceleration to direct COM acceleration",
        y=0.995,
        fontsize=15,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(OUT / "gain_only_model_vs_measured.png", dpi=200)
    plt.close(fig)

    pd.DataFrame(summaries).to_csv(
        OUT / "identified_parameters_and_fit.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(frequency_rows).to_csv(
        OUT / "frequency_response_values.csv", index=False, encoding="utf-8-sig"
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
