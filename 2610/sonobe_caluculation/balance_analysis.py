"""Analysis of MedicalForcePlate's processed Sonobe strategy-mode signals.

Reference: Sonobe & Mitsutake (2025), doi:10.3389/fneur.2025.1500308,
sections 2.2, 2.3, 2.5. No raw-IMU attitude estimation is performed here.
"""
from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REFERENCE = 'https://doi.org/10.3389/fneur.2025.1500308'
MODE_COLUMNS = ('xix1', 'xix2', 'xiy1', 'xiy2')
AXES = {'AP': ('xix1', 'xix2', 'axhd'), 'ML': ('xiy1', 'xiy2', 'ayhd')}


@dataclass
class Recording:
    path: Path
    time: np.ndarray
    signals: dict[str, np.ndarray]
    metadata: dict[str, str]
    exported_indices: dict[str, float]
    encoding: str
    fs: float


def read_recording(path: str | Path) -> Recording:
    path = Path(path)
    raw = path.read_bytes()
    for encoding in ('utf-8-sig', 'cp932'):
        try:
            text = raw.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise ValueError('CSVの文字コードを読めません。UTF-8またはShift-JISで保存してください。')
    rows = list(csv.reader(io.StringIO(text)))
    # MedicalForcePlate exports the same header twice (statistics, then samples).
    headers = [i for i, row in enumerate(rows) if set(MODE_COLUMNS).issubset(row)]
    if not headers:
        raise ValueError('処理済みモード列 xix1, xix2, xiy1, xiy2 が必要です。'
                         '生の床反力・IMUだけのCSVは、この版では解析できません。')
    index = headers[-1]
    header = [s.strip() for s in rows[index]]
    if len(set(header)) != len(header):
        raise ValueError('時系列ヘッダーに重複した列名があります。')
    numeric = []
    for line, row in enumerate(rows[index + 1:], start=index + 2):
        if not row or not any(s.strip() for s in row):
            continue
        if len(row) != len(header):
            raise ValueError(f'{line}行目: 列数が一致しません。')
        # Only consume required numeric columns; unused placeholders may be blank.
        needed = [0] + [header.index(c) for c in (*MODE_COLUMNS, 'axhd', 'ayhd') if c in header]
        try:
            values = [float(row[k]) for k in needed]
        except ValueError as exc:
            raise ValueError(f'{line}行目: 時刻または必要な加速度が数値ではありません。') from exc
        if not np.isfinite(values).all():
            raise ValueError(f'{line}行目: 欠損値または無限大を含みます。補間せず停止しました。')
        numeric.append(values)
    if len(numeric) < 3:
        raise ValueError('時系列データが3点未満です。')
    data = np.asarray(numeric)
    time = data[:, 0]
    dt = np.diff(time)
    if np.any(dt <= 0):
        raise ValueError('時刻は重複のない昇順である必要があります。')
    interval = float(np.median(dt))
    if not np.allclose(dt, interval, rtol=0.01, atol=1e-7):
        raise ValueError('時刻間隔が不均一です。欠落・同期を確認してください。')
    names = [c for c in (*MODE_COLUMNS, 'axhd', 'ayhd') if c in header]
    metadata = {r[0]: r[1] for r in rows[:headers[0]] if len(r) == 2}
    if '記録周波数(Hz)' in metadata:
        declared = float(metadata['記録周波数(Hz)'])
        if not np.isclose(declared, 1 / interval, rtol=0.01):
            raise ValueError('記録周波数と時刻列が一致しません。')
    # Validate units when present; compact custom CSVs must use documented mm/s².
    unit_rows = [r for r in rows[:index] if r and r[0] == '単位' and len(r) == len(header)]
    if unit_rows:
        for c in names:
            if unit_rows[-1][header.index(c)] not in ('mm/s/s', 'mm/s^2', 'mm/s²'):
                raise ValueError(f'{c}の単位はmm/s²である必要があります。')
    indices = {}
    for i, row in enumerate(rows[:index]):
        if 'RMSXI1ap' in row and i + 2 < len(rows):
            try:
                indices = {key: float(value) for key, value in zip(row[1:], rows[i+2][1:])}
            except ValueError:
                indices = {}
    return Recording(path, time, dict(zip(names, data[:, 1:].T)), metadata,
                     indices, encoding, 1 / interval)


def mode_metrics(ankle: np.ndarray, hip: np.ndarray) -> dict:
    ankle, hip = np.asarray(ankle, float), np.asarray(hip, float)
    if ankle.ndim != 1 or ankle.shape != hip.shape or len(ankle) < 3:
        raise ValueError('同じ長さの1次元波形が3点以上必要です。')
    if not np.isfinite(ankle).all() or not np.isfinite(hip).all():
        raise ValueError('波形に欠損があります。')
    a, h = ankle - ankle.mean(), hip - hip.mean()
    ra, rh = float(np.sqrt(np.mean(a*a))), float(np.sqrt(np.mean(h*h)))
    total = ra + rh
    # Undefined ratios remain None, never silently become zero or 50%.
    return dict(ankle_rms=ra, hip_rms=rh,
                ankle_percent=100 * ra / total if total > 1e-12 else None,
                hip_percent=100 * rh / total if total > 1e-12 else None,
                rxi=rh / ra if ra > 1e-12 else None,
                cxi=float(np.clip(np.mean(a*h)/(ra*rh), -1, 1)) if min(ra, rh) > 1e-12 else None)


def analyze(recording: Recording, start: float = 5, end: float = 35) -> dict:
    t = recording.time
    if not np.isfinite([start, end]).all() or start >= end:
        raise ValueError('開始・終了は有限の秒数で、開始 < 終了にしてください。')
    if start < t[0] - 1e-7 or end > t[-1] + 1/recording.fs + 1e-7:
        raise ValueError(f'解析区間が記録範囲 {t[0]:g}〜{t[-1]+1/recording.fs:g} 秒を超えています。')
    mask = (t >= start - 1e-9) & (t < end - 1e-9)
    if mask.sum() < 3:
        raise ValueError('解析区間に3点以上のデータが必要です。')
    result = dict(source=recording.path.name, reference=REFERENCE, start_s=start, end_s=end,
                  interval='start inclusive, end exclusive', samples=int(mask.sum()),
                  sampling_hz=recording.fs, input='exported_processed_modes',
                  percentage_definition='100 * RMS(mode) / (RMS(ankle) + RMS(hip)); demeaned',
                  additional_filter=False, axes={}, warnings=[])
    for axis, (ca, ch, head) in AXES.items():
        a, h = recording.signals[ca][mask], recording.signals[ch][mask]
        metrics = mode_metrics(a, h)
        if metrics['ankle_percent'] is None:
            result['warnings'].append(f'{axis}: 両モードの変動がゼロのため割合は未定義です。')
        if metrics['rxi'] is None:
            result['warnings'].append(f'{axis}: 足関節RMSがゼロのためRXIは未定義です。')
        if head in recording.signals:
            hd = recording.signals[head][mask]
            residual = (a-a.mean()) + (h-h.mean()) - (hd-hd.mean())
            err = float(np.sqrt(np.mean(residual**2)))
            hrms = float(np.std(hd))
            metrics.update(head_rms=hrms, reconstruction_rmse=err,
                           reconstruction_relative_percent=100*err/hrms if hrms > 1e-12 else None)
            # Engineering consistency check, not a clinical threshold.
            if err > max(1e-5, 0.01 * hrms):
                result['warnings'].append(
                    f'{axis}: モード和と頭部加速度にずれがあります '
                    f'(平均除去後RMSE {err:.3f} mm/s²)。モード割合は出力列から計算していますが、'
                    '頭部加速度への寄与の解釈には機器側の計算定義確認が必要です。')
        else:
            result['warnings'].append(f'{axis}: 頭部加速度列がないためモード和の整合性は未確認です。')
        # Comparison only; an exported summary may use a different analysis interval.
        suffix = axis.lower()
        metrics['exported_summary'] = {k: recording.exported_indices.get(col + suffix)
                                      for k, col in [('ankle_rms','RMSXI1'), ('hip_rms','RMSXI2'),
                                                     ('rxi','RXI'), ('cxi','CXI')]}
        result['axes'][axis] = metrics
    return result


def make_figure(recording: Recording, result: dict):
    from matplotlib.figure import Figure
    from matplotlib import rcParams
    rcParams['font.family'] = ['Yu Gothic', 'Meiryo', 'DejaVu Sans']
    figure = Figure(figsize=(12, 8.4), dpi=110, facecolor='#f5f7fb', layout='constrained')
    grid = figure.add_gridspec(3, 2, height_ratios=[1.1, 1.8, 0.8])
    figure.suptitle(f"足関節・股関節モード  |  {recording.path.name}\n"
                   f"{result['start_s']:g} ≤ t < {result['end_s']:g} 秒 · {result['samples']:,} samples",
                   fontsize=16, fontweight='bold', color='#172c45')
    footer = '割合はRMS合計による正規化（時間・トルク・分散の寄与率ではありません）。'
    if result['warnings']:
        footer += '\n注意あり：モード和と頭部加速度の不一致・未定義等。詳細は画面 / notes.txt を確認。'
    figure.supxlabel(footer, fontsize=9, color='#844710')
    mask = (recording.time >= result['start_s'] - 1e-9) & (recording.time < result['end_s'] - 1e-9)
    colors = ['#167d9a', '#e08a3b']
    for i, (axis, (ca, ch, head)) in enumerate(AXES.items()):
        m = result['axes'][axis]
        bar = figure.add_subplot(grid[0, i])
        bar.set_title('前後方向 (AP)' if axis == 'AP' else '左右方向 (ML)', loc='left', fontweight='bold')
        a, h = m['ankle_percent'], m['hip_percent']
        if a is not None:
            bar.barh(0, a, color=colors[0], height=.5)
            bar.barh(0, h, left=a, color=colors[1], height=.5)
            bar.text(.01, .05, f'足関節 {a:.1f}%    股関節 {h:.1f}%', transform=bar.transAxes, fontsize=13)
        else:
            bar.text(50, 0, '変動ゼロ：割合は未定義', ha='center')
        bar.set(xlim=(0, 100), ylim=(-.6, .6), yticks=[], xlabel='RMS合計に対する割合 (%)')
        for spine in bar.spines.values():
            spine.set_visible(False)
        trace = figure.add_subplot(grid[1, i])
        for column, label, color in [(ca, '足関節', colors[0]), (ch, '股関節', colors[1]),
                                      (head, '頭部', '#697a8e')]:
            if column in recording.signals:
                signal = recording.signals[column][mask]
                trace.plot(recording.time[mask], signal-signal.mean(), label=label,
                           color=color, linewidth=.85, alpha=.85)
        trace.set(xlabel='時刻 (s)', ylabel='平均除去後の加速度 (mm/s²)')
        trace.grid(alpha=.17)
        trace.legend(loc='upper right', ncol=3, fontsize=9)
        summary = figure.add_subplot(grid[2, i])
        summary.axis('off')
        fmt = lambda v: '未定義' if v is None else f'{v:.3f}'
        summary.text(0, .95, f"足関節RMS   {m['ankle_rms']:.3f} mm/s²\n"
                     f"股関節RMS   {m['hip_rms']:.3f} mm/s²\n"
                     f"RXI (股 / 足)   {fmt(m['rxi'])}     CXI   {fmt(m['cxi'])}",
                     va='top', fontsize=11, linespacing=1.6)
    return figure


def export_results(recording: Recording, result: dict, directory: str | Path) -> Path:
    from datetime import datetime
    # A new directory per export keeps prior results intact.
    target = Path(directory) / (recording.path.stem + '_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    target.mkdir(parents=True, exist_ok=False)
    (target/'result.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    with (target/'summary.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['source', 'start_s', 'end_s', 'samples', 'axis', 'ankle_rms_mm_s2',
                         'hip_rms_mm_s2', 'ankle_rms_share_percent', 'hip_rms_share_percent', 'RXI', 'CXI'])
        for axis, m in result['axes'].items():
            writer.writerow([result['source'], result['start_s'], result['end_s'], result['samples'], axis,
                             *[m[k] for k in ('ankle_rms','hip_rms','ankle_percent','hip_percent','rxi','cxi')]])
    mask = (recording.time >= result['start_s']-1e-9) & (recording.time < result['end_s']-1e-9)
    with (target/'waveforms.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        names = list(recording.signals)
        writer.writerow(['time_s']+[c+'_demeaned_mm_s2' for c in names])
        signals = [recording.signals[c][mask] - recording.signals[c][mask].mean() for c in names]
        writer.writerows(zip(recording.time[mask], *signals))
    fig = make_figure(recording, result)
    fig.savefig(target/'result.png', dpi=160)
    notes = ('割合 = 各モードの平均除去後RMS / 両モードRMSの合計 × 100。\n'
             '百分率は表示用の正規化であり、論文のRXIそのものではありません。\n'
             '時間割合・関節角度割合・トルク割合・分散寄与率ではありません。\n'
             'CSVの処理済みモードを使用。追加フィルタなし。生IMUの姿勢推定は実行していません。\n'
             'MLの股関節モードは論文の腰部関節モデルに対応します。\n'
             '整合性検査の1%閾値は実装上の目安であり、臨床判定基準ではありません。\n'
             f'参考文献: {REFERENCE}\n\n' + '\n'.join(result['warnings']))
    (target/'notes.txt').write_text(notes, encoding='utf-8-sig')
    return target
