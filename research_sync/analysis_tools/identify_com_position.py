"""Fit the visual path from logged background velocity to estimated body COM."""
import argparse
import json
from pathlib import Path

import numpy as np

from analyze_periodic_response import harmonic, metrics, period, table, write
from fit_visual_path import G, base_response, profile


def cycle_fft(t, signal, start, duration, frequencies):
    """Return complex peak amplitudes at the exact periodic stimulus bins."""
    count = round(100 * duration)
    grid = start + np.arange(count) / 100
    valid = np.isfinite(t) & np.isfinite(signal)
    if grid[0] < t[valid][0] or grid[-1] > t[valid][-1]:
        raise ValueError('FFT grid outside valid observations')
    values = np.interp(grid, t[valid], signal[valid])
    bins = np.rint(frequencies * duration).astype(int)
    if not np.allclose(bins, frequencies * duration):
        raise ValueError('Stimulus frequency outside the periodic FFT bins')
    return 2 * np.fft.rfft(values)[bins] / count


def analyze(folder, fft=False, cycle_data=None, time_data=None, reaction_sign=-1):
    stem = folder.name
    meta = json.loads((folder / f'{stem}_metadata.json').read_text(encoding='utf-8-sig'))
    stim = meta['stimulus']
    f = np.array(stim['frequencies_hz'])
    if meta['termination'] != 'complete' or np.unique(f[(f > 0) & (f < .5)]).size < 2:
        raise ValueError(f'{stem}: expected a complete trial with two or more stimulus frequencies below 0.5 Hz')
    start, stop = stim['analysis_interval_s']
    cycle_s = period(f)
    cycles = round((stop - start) / cycle_s)
    assert np.isclose(stop - start, cycles * cycle_s)
    q, _ = table(folder / f'{stem}_quest.csv')
    b, _ = table(folder / f'{stem}_combalance.csv')
    qt, qe, ft = q['pc_mapped_time_s'], q['stimulus_elapsed_s'], b['pc_estimated_sample_time_s']
    for t in (qt, qe, ft):
        if not np.isfinite(t).all() or not np.all(np.diff(t) > 0):
            raise ValueError(f'{stem}: invalid timestamps')
    if (qe[0] > start or qe[-1] < stop or
            np.any(np.interp([start, stop], qe, qt) < ft[0]) or
            np.any(np.interp([start, stop], qe, qt) > ft[-1]) or
            max(np.diff(qt).max(), np.diff(ft).max()) >= .1):
        raise ValueError(f'{stem}: insufficient contiguous coverage')
    elapsed = np.interp(ft, qt, qe, left=np.nan, right=np.nan)
    starts = start + cycle_s * np.arange(cycles)
    force_masks = [(elapsed >= a) & (elapsed < a + cycle_s) for a in starts]
    u = []
    for a in start + cycle_s * np.arange(cycles):
        mask = (qe >= a) & (qe < a + cycle_s)
        u.append(cycle_fft(qe, q['ap_velocity_m_s'], a, cycle_s, f) if fft
                 else harmonic(qe[mask], q['ap_velocity_m_s'][mask], f))
    u = np.array(u)
    keep = f < .5
    mass = meta['weight_estimate']['body_mass_kg']
    height = meta['participant']['height_cm'] / 100
    mb = .978 * mass
    # OEM AP convention from standing_analysis: plate force Fy, AP moment Mx.
    # Use reaction on the body in the user equation. COP uses the recorded
    # moment origin and instantaneous vertical load (no contact-plane shift).
    rx = b['sensor_Fy_N']
    vertical = b['sensor_Fz_N']
    valid = ((elapsed >= start - .1) & (elapsed < stop + .1) if fft
             else np.logical_or.reduce(force_masks))
    if np.any(vertical[valid] < 20):
        raise ValueError(f'{stem}: insufficient standing load')
    xp = b['sensor_Mx_Nm'][valid] / vertical[valid]
    fits, points = [], []
    for anthro, (j, l) in {'male': (.0473, .513), 'female': (.0481, .504)}.items():
        lb = l * height
        jb = j * mass * height**2
        c = jb / (mb * lb) + lb + .039 * height
        if time_data is not None:
            loaded = np.isfinite(elapsed) & (vertical >= 20)
            position = mass / mb * b['sensor_Mx_Nm'][loaded] / vertical[loaded]
            position += reaction_sign * c / (mb * G) * rx[loaded]
            time_data[anthro] = (qe, q['ap_displacement_m'], q['ap_velocity_m_s'],
                                 elapsed[loaded], position)
        xb = np.full(len(ft), np.nan)
        xb[valid] = mass / mb * xp + reaction_sign * c / (mb * G) * rx[valid]
        y = np.array([cycle_fft(elapsed, xb, a, cycle_s, f) for a in starts] if fft
                     else [harmonic(elapsed[mask], xb[mask], f) for mask in force_masks])
        response, _, snr, _, coherence, *_ = metrics(u, y)
        cross = np.mean(np.conj(u) * y, axis=0)
        input_power = np.mean(abs(u)**2, axis=0)
        output_power = np.mean(abs(y)**2, axis=0)
        np.testing.assert_allclose(response, cross / input_power)
        s = 2j * np.pi * f
        base = base_response(f, mass, height, anthro, 'full') / ((mb / mass) * s**2)
        if cycle_data is not None:
            cycle_data[anthro] = (f[keep], u[:, keep], y[:, keep], base[keep])
        for mode in ('gain_only', 'complex'):
            summary, _, _, pred = profile(f[keep], response[keep], base[keep], mode)
            sensitivity = []
            for excluded in range(cycles):
                included = np.arange(cycles) != excluded
                h = np.sum(np.conj(u[included]) * y[included], axis=0) / np.sum(abs(u[included])**2, axis=0)
                alternate, *_ = profile(f[keep], h[keep], base[keep], mode)
                sensitivity.append(alternate['Tv_s'])
            fits.append(dict(recording=stem, condition=stim['condition'],
                             anthropometry=anthro, fit=mode,
                             leave_one_cycle_out_Tv_s=';'.join(f'{v:.6g}' for v in sensitivity),
                             median_coherence=float(np.median(coherence[keep])), **summary))
            for i, observed, predicted in zip(np.flatnonzero(keep), response[keep], pred):
                points.append(dict(recording=stem, anthropometry=anthro, fit=mode,
                                   frequency_hz=f[i], observed_gain=abs(observed),
                                   observed_phase_deg=np.angle(observed, deg=True),
                                   predicted_gain=abs(predicted),
                                   predicted_phase_deg=np.angle(predicted, deg=True),
                                   coherence=coherence[i], snr_db=snr[i]))
                points[-1].update(cross_real=cross[i].real, cross_imag=cross[i].imag,
                                  input_power=input_power[i], output_power=output_power[i])
    return fits, points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trials', nargs='+', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--fft', action='store_true', help='Use rectangular full-cycle FFT coefficients')
    args = parser.parse_args()
    fits, points = [], []
    for folder in args.trials:
        rows, frf = analyze(folder, fft=args.fft)
        fits.extend(rows)
        points.extend(frf)
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / 'fits.csv', fits)
    write(args.output / 'frequency_response.csv', points)


if __name__ == '__main__':
    main()
