from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


SOURCE_ROOT = Path(r"C:\Users\seguchi\Desktop\1027")
OUTPUT_PATH = Path(r"C:\Users\seguchi\gitseguchi\All\2610\outputs\fft_1027_20261101\fft_data.json")


def analyze_file(path: Path) -> dict:
    timestamps_ms: list[float] = []
    accel_z: list[float] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 5 or row[0] != "ags":
                continue
            timestamps_ms.append(float(row[1]))
            accel_z.append(float(row[4]))

    timestamps = np.asarray(timestamps_ms, dtype=float)
    values = np.asarray(accel_z, dtype=float)
    relative_seconds = (timestamps - timestamps[0]) / 1000.0
    mask = (relative_seconds >= 5.0) & (relative_seconds < 30.0)
    segment_time = relative_seconds[mask]
    segment = values[mask]

    dt_seconds = float(np.median(np.diff(timestamps)) / 1000.0)
    sample_rate_hz = 1.0 / dt_seconds
    centered = segment - segment.mean()
    frequencies = np.fft.rfftfreq(centered.size, d=dt_seconds)
    amplitudes = 2.0 * np.abs(np.fft.rfft(centered)) / centered.size
    if amplitudes.size:
        amplitudes[0] /= 2.0
    band = frequencies <= 1.0 + 1e-12

    return {
        "condition": path.parent.name.removeprefix("20261027-"),
        "source_file": str(path),
        "total_samples": int(values.size),
        "source_duration_seconds": float(relative_seconds[-1]),
        "sample_interval_seconds": dt_seconds,
        "sample_rate_hz": sample_rate_hz,
        "window_start_seconds": float(segment_time[0]),
        "window_end_seconds": float(segment_time[-1]),
        "window_samples": int(segment.size),
        "frequency_resolution_hz": float(sample_rate_hz / segment.size),
        "frequencies_hz": frequencies[band].round(12).tolist(),
        "amplitudes": amplitudes[band].tolist(),
    }


def main() -> None:
    files = sorted(
        path
        for path in SOURCE_ROOT.rglob("*.csv")
        if "なんでもない" not in str(path.parent)
    )
    results = [analyze_file(path) for path in files]
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(
        json.dumps({"method": "mean-removed one-sided FFT amplitude", "results": results}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    for item in results:
        peak_index = int(np.argmax(item["amplitudes"][1:]) + 1)
        print(
            f"{item['condition']}: total={item['total_samples']}, window={item['window_samples']}, "
            f"fs={item['sample_rate_hz']:.6f} Hz, df={item['frequency_resolution_hz']:.6f} Hz, "
            f"peak(>0)={item['frequencies_hz'][peak_index]:.2f} Hz / {item['amplitudes'][peak_index]:.6f}"
        )


if __name__ == "__main__":
    main()
