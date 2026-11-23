"""Compare body-sway magnitude across three 30-s fixation conditions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import detrend, windows


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "fixation_condition_sway_comparison"
CONDITION_LABELS = {
    "none": "No fixation point",
    "dot": "Red fixed point",
    "spinner": "Rotating loading ring",
}
COLORS = {
    "none": "#4C78A8",
    "dot": "#E45756",
    "spinner": "#59A14F",
}
STIMULUS_FREQUENCY_HZ = 0.20
ANALYSIS_START_S = 3.0
ANALYSIS_END_S = 27.0


def sinusoid_amplitude(t: np.ndarray, y: np.ndarray, frequency: float) -> float:
    """Least-squares sinusoidal amplitude after removing offset and drift."""
    t = np.asarray(t, float) - float(t[0])
    y = np.asarray(y, float)
    nuisance = np.column_stack((np.ones_like(t), t))
    y = y - nuisance @ np.linalg.lstsq(nuisance, y, rcond=None)[0]
    omega = 2 * np.pi * frequency
    basis = np.column_stack((np.cos(omega * t), np.sin(omega * t)))
    cosine, sine = np.linalg.lstsq(basis, y, rcond=None)[0]
    return float(np.hypot(cosine, sine))


def amplitude_spectrum(t: np.ndarray, y: np.ndarray):
    sample_interval = float(np.median(np.diff(t)))
    sample_rate = 1.0 / sample_interval
    y = detrend(np.asarray(y, float), type="linear")
    taper = windows.hann(len(y), sym=False)
    fft_length = max(16384, int(2 ** np.ceil(np.log2(len(y)))))
    spectrum = np.fft.rfft(y * taper, n=fft_length)
    frequencies = np.fft.rfftfreq(fft_length, sample_interval)
    amplitude = 2 * np.abs(spectrum) / np.sum(taper)
    return frequencies, amplitude


def load_condition(folder: Path):
    stem = folder.name
    quest = pd.read_csv(folder / f"{stem}_quest.csv")
    force = pd.read_csv(folder / f"{stem}_combalance.csv")
    condition = str(quest["visual_condition"].iloc[0])
    total_mass = float(quest["body_mass_kg"].iloc[0])
    height = float(quest["height_cm"].iloc[0]) / 100

    # Male anthropometric coefficients supplied in the visual-balance note.
    body_mass = 0.978 * total_mass
    body_inertia = 0.0425 * total_mass * height**2
    com_height = 0.531 * height
    ankle_height = 0.038 * height
    force_coefficient = (
        body_inertia / (body_mass * com_height) + com_height + ankle_height
    )

    quest_time = quest["pc_mapped_time_s"].to_numpy(float)
    quest_elapsed = quest["stimulus_elapsed_s"].to_numpy(float)
    start_pc, end_pc = np.interp(
        [ANALYSIS_START_S, ANALYSIS_END_S], quest_elapsed, quest_time
    )
    force_time = force["pc_estimated_sample_time_s"].to_numpy(float)
    keep = (force_time >= start_pc) & (force_time < end_pc)
    time = force_time[keep] - start_pc

    # Participant AP axis after the 90-degree CCW turn.
    horizontal_force = force.loc[keep, "sensor_Fy_N"].to_numpy(float)
    sagittal_moment = -force.loc[keep, "sensor_Mx_Nm"].to_numpy(float)
    com_acceleration = horizontal_force / total_mass
    com_displacement = (
        sagittal_moment - force_coefficient * horizontal_force
    ) / body_mass

    displacement_detrended = detrend(com_displacement, type="linear")
    acceleration_detrended = detrend(com_acceleration, type="linear")
    return {
        "condition": condition,
        "label": CONDITION_LABELS[condition],
        "color": COLORS[condition],
        "time": time,
        "com_displacement": displacement_detrended,
        "com_acceleration": acceleration_detrended,
        "displacement_rms_m": float(np.sqrt(np.mean(displacement_detrended**2))),
        "displacement_amplitude_0p2_m": sinusoid_amplitude(
            time, com_displacement, STIMULUS_FREQUENCY_HZ
        ),
        "acceleration_rms_m_s2": float(np.sqrt(np.mean(acceleration_detrended**2))),
        "acceleration_amplitude_0p2_m_s2": sinusoid_amplitude(
            time, com_acceleration, STIMULUS_FREQUENCY_HZ
        ),
        "samples": int(keep.sum()),
    }


def main():
    OUT.mkdir(exist_ok=True)
    folders = sorted(ROOT.glob("20261121_otsubo_30s_*_001"))
    records = [load_condition(folder) for folder in folders]
    order = ["none", "dot", "spinner"]
    records.sort(key=lambda record: order.index(record["condition"]))

    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    displacement_time_ax, displacement_fft_ax, displacement_bar_ax = axes[0]
    acceleration_time_ax, acceleration_fft_ax, acceleration_bar_ax = axes[1]

    for record in records:
        time = record["time"]
        color = record["color"]
        label = record["label"]
        displacement_time_ax.plot(
            time, 1000 * record["com_displacement"], color=color, label=label, linewidth=1.2
        )
        acceleration_time_ax.plot(
            time, record["com_acceleration"], color=color, label=label, linewidth=1.2
        )
        frequency, displacement_spectrum = amplitude_spectrum(
            time, record["com_displacement"]
        )
        _, acceleration_spectrum = amplitude_spectrum(
            time, record["com_acceleration"]
        )
        show = frequency <= 1.0
        displacement_fft_ax.plot(
            frequency[show], 1000 * displacement_spectrum[show], color=color, label=label
        )
        acceleration_fft_ax.plot(
            frequency[show], acceleration_spectrum[show], color=color, label=label
        )

    for ax in (displacement_fft_ax, acceleration_fft_ax):
        ax.axvline(
            STIMULUS_FREQUENCY_HZ,
            color="black",
            linestyle="--",
            linewidth=1,
            label="Stimulus 0.20 Hz",
        )
        ax.set_xlim(0, 1.0)

    labels = [record["label"] for record in records]
    colors = [record["color"] for record in records]
    x = np.arange(len(records))
    width = 0.36
    displacement_bar_ax.bar(
        x - width / 2,
        [1000 * record["displacement_amplitude_0p2_m"] for record in records],
        width,
        color=colors,
        alpha=0.95,
        label="0.20-Hz amplitude",
    )
    displacement_bar_ax.bar(
        x + width / 2,
        [1000 * record["displacement_rms_m"] for record in records],
        width,
        color=colors,
        alpha=0.45,
        hatch="//",
        label="RMS (3-27 s)",
    )
    acceleration_bar_ax.bar(
        x - width / 2,
        [record["acceleration_amplitude_0p2_m_s2"] for record in records],
        width,
        color=colors,
        alpha=0.95,
        label="0.20-Hz amplitude",
    )
    acceleration_bar_ax.bar(
        x + width / 2,
        [record["acceleration_rms_m_s2"] for record in records],
        width,
        color=colors,
        alpha=0.45,
        hatch="//",
        label="RMS (3-27 s)",
    )

    displacement_time_ax.set_title("Estimated COM displacement")
    displacement_time_ax.set_xlabel("Time in steady interval [s]")
    displacement_time_ax.set_ylabel("Detrended displacement [mm]")
    acceleration_time_ax.set_title("Direct COM acceleration from horizontal force")
    acceleration_time_ax.set_xlabel("Time in steady interval [s]")
    acceleration_time_ax.set_ylabel("Detrended acceleration [m/s²]")
    displacement_fft_ax.set_title("COM displacement spectrum")
    displacement_fft_ax.set_xlabel("Frequency [Hz]")
    displacement_fft_ax.set_ylabel("Amplitude [mm]")
    acceleration_fft_ax.set_title("COM acceleration spectrum")
    acceleration_fft_ax.set_xlabel("Frequency [Hz]")
    acceleration_fft_ax.set_ylabel("Amplitude [m/s²]")
    displacement_bar_ax.set_title("COM displacement magnitude")
    displacement_bar_ax.set_ylabel("Magnitude [mm]")
    acceleration_bar_ax.set_title("COM acceleration magnitude")
    acceleration_bar_ax.set_ylabel("Magnitude [m/s²]")

    for ax in axes.ravel():
        ax.grid(True, alpha=0.25)
    for ax in (displacement_bar_ax, acceleration_bar_ax):
        ax.set_xticks(x, ["None", "Red dot", "Loading ring"], rotation=12)
        ax.legend(frameon=False, fontsize=8)
    displacement_time_ax.legend(frameon=False, fontsize=8)
    acceleration_time_ax.legend(frameon=False, fontsize=8)
    displacement_fft_ax.legend(frameon=False, fontsize=8)
    acceleration_fft_ax.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Otsubo: body-sway magnitude across fixation-display conditions (steady 3-27 s)",
        fontsize=15,
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "fixation_condition_sway_comparison.png", dpi=200)
    plt.close(fig)

    summary = pd.DataFrame(
        [
            {
                "condition": record["condition"],
                "condition_label": record["label"],
                "analysis_interval_stimulus_s": "3-27",
                "samples": record["samples"],
                "COM_displacement_RMS_mm": 1000 * record["displacement_rms_m"],
                "COM_displacement_0p20Hz_amplitude_mm": 1000
                * record["displacement_amplitude_0p2_m"],
                "COM_acceleration_RMS_m_s2": record["acceleration_rms_m_s2"],
                "COM_acceleration_0p20Hz_amplitude_m_s2": record[
                    "acceleration_amplitude_0p2_m_s2"
                ],
            }
            for record in records
        ]
    )
    baseline = summary.iloc[0]
    for metric in (
        "COM_displacement_RMS_mm",
        "COM_displacement_0p20Hz_amplitude_mm",
        "COM_acceleration_RMS_m_s2",
        "COM_acceleration_0p20Hz_amplitude_m_s2",
    ):
        summary[f"{metric}_relative_to_no_fixation_percent"] = (
            100 * summary[metric] / baseline[metric]
        )
    summary.to_csv(OUT / "fixation_condition_sway_summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
