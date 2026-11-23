"""Plot FFT spectra of the seven external (head) IMU channels."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "head_imu_fft"
CHANNELS = [f"external_imu_count_{i}" for i in range(7)]
INPUT_FREQUENCIES = np.arange(0.05, 1.5000001, 0.05)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    folders = sorted(ROOT.glob("20261121_*_90s_multicosine_*"))
    summary_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, folder in zip(axes, folders):
        stem = folder.name
        data = pd.read_csv(folder / f"{stem}_combalance.csv")
        fs = float(data["sampling_hz"].median())
        # The force/IMU recording starts 10 s before the visual stimulus. The
        # pure multicosine interval is stimulus elapsed 10-90 s, hence force
        # recording elapsed 20-100 s.
        elapsed = data["sample_elapsed_s"].to_numpy(float)
        keep = (elapsed >= 20.0) & (elapsed < 100.0)
        values = data.loc[keep, CHANNELS].to_numpy(float)
        n = len(values)
        frequencies = np.fft.rfftfreq(n, 1.0 / fs)
        demeaned = values - np.mean(values, axis=0, keepdims=True)
        amplitudes = 2.0 * np.abs(np.fft.rfft(demeaned, axis=0)) / n

        band = (frequencies > 0) & (frequencies <= 2.0)
        for input_frequency in INPUT_FREQUENCIES:
            ax.axvline(input_frequency, color="0.82", linewidth=0.55, zorder=0)
        for index, channel in enumerate(CHANNELS):
            ax.plot(
                frequencies[band], amplitudes[band, index],
                linewidth=1.0, label=channel.replace("external_imu_count_", "ch"),
            )
            non_dc = (frequencies >= 0.01) & (frequencies <= 2.0)
            peak_index = np.argmax(amplitudes[non_dc, index])
            peak_frequency = frequencies[non_dc][peak_index]
            peak_amplitude = amplitudes[non_dc, index][peak_index]
            summary_rows.append(
                {
                    "trial": stem,
                    "channel": channel,
                    "samples": n,
                    "sampling_hz": fs,
                    "minimum_count": np.min(values[:, index]),
                    "maximum_count": np.max(values[:, index]),
                    "nonzero_samples": int(np.count_nonzero(values[:, index])),
                    "peak_frequency_hz_0_to_2": peak_frequency,
                    "peak_amplitude_count": peak_amplitude,
                }
            )

        subject = str(data["subject"].iloc[0])
        trial = int(data["trial_number"].iloc[0])
        ax.set_title(f"{subject}, trial {trial}")
        ax.set_xlim(0, 2.0)
        ax.grid(axis="y", alpha=0.25)
        ax.set_ylabel("FFT amplitude [count]")

    for ax in axes[-2:]:
        ax.set_xlabel("Frequency [Hz]")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=7, frameon=False)
    fig.suptitle(
        "Head IMU FFT (gray vertical lines: visual-input frequencies 0.05-1.50 Hz)",
        y=0.97,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT / "head_imu_fft_all_trials.png", dpi=180)
    plt.close(fig)

    pd.DataFrame(summary_rows).to_csv(
        OUT / "head_imu_fft_peaks.csv", index=False, encoding="utf-8-sig"
    )
    print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
