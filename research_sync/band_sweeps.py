"""Six 100-second AP sweeps; constant velocity envelope in the analysis core."""
import argparse
import json
from pathlib import Path

EDGES = tuple(.05 * 10 ** (i / 3) for i in range(4))
EDGES = (*EDGES[:3], .5)
PROFILES = {}
LABELS = {}
FILENAMES = {}
for band in range(3):
    lo, hi = EDGES[band:band + 2]
    for direction in ('up', 'down'):
        condition = f'study_sweep100_band{band+1}_{direction}'
        start, end = (lo, hi) if direction == 'up' else (hi, lo)
        PROFILES[condition] = dict(
            schema=8, waveform='ap_velocity_log_band_sweep', condition=condition,
            enabled=True, duration_s=100, onset_s=10, sweep_s=85, offset_s=5,
            analysis_interval_s=[10, 95], band=band+1, direction=direction,
            start_frequency_hz=start, end_frequency_hz=end,
            velocity_amplitude_m_s=.1, sample_step_s=.001,
            initial_phase_rad=0, initial_position_m=0,
            hold_final_position=True, fixation=None,
            envelope='quintic smoothstep velocity envelope; fixed frequency in onset/offset',
            phase='continuous phase; exponential frequency in 10 <= t <= 95 seconds',
            position='trapezoidal integral of velocity at 1 ms; no return to origin',
            analysis_note='Onset/offset excluded; 10-second onset does not establish stationarity',
        )
        arrow = '上昇' if direction == 'up' else '下降'
        LABELS[f'100秒 分割sweep{band+1} {arrow} {start:.3f}→{end:.3f} Hz'] = condition
        FILENAMES[condition] = f'100s_sweep_band{band+1}_{direction}'

def profile(condition):
    # Return a fresh object, including its nested analysis interval.
    return json.loads(json.dumps(PROFILES[condition]))

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--export',type=Path,default=Path(__file__).resolve().parent/'study_profiles')
    args=parser.parse_args();args.export.mkdir(parents=True,exist_ok=True)
    for key in PROFILES:
        (args.export/f'{key}.json').write_text(json.dumps(profile(key),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(f'Exported six profiles to {args.export}')
