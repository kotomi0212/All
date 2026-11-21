from pathlib import Path
import csv
import html
import shutil

SOURCE = Path(r"C:\Users\seguchi\B4春\261023\head_imu_20261026")
OUT = Path(__file__).resolve().parent / "head_imu_20261026_AP_only"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def display_name(stem):
    return stem.removeprefix("261022_").replace("_seguchi", " / ")


def category(stem):
    key = stem.lower()
    if "statics" in key:
        return "静止"
    if "hip" in key:
        return "腰部条件"
    if "leg" in key:
        return "下肢条件"
    if "trans" in key:
        return "移動条件"
    return "その他"


def polyline(values, x0, y0, width, height, ymin=None, ymax=None):
    if ymin is None:
        ymin = min(values)
    if ymax is None:
        ymax = max(values)
    span = ymax - ymin or 1
    step = max(1, len(values) // 1800)
    points = []
    for i in range(0, len(values), step):
        x = x0 + width * i / max(1, len(values) - 1)
        y = y0 + height * (ymax - values[i]) / span
        points.append(f"{x:.1f},{y:.1f}")
    return " ".join(points)


def record_svg(path, title, t, ap, freq, amp, rms, full_rms, peak_hz):
    shown = [v for f, v in zip(freq, amp) if f <= 6]
    wave_pts = polyline(ap, 75, 105, 1065, 245)
    fft_pts = polyline(shown, 75, 465, 1065, 245, 0, max(shown) * 1.08)
    shade_x, shade_w = 75 + 1065 * 5 / 60, 1065 * 30 / 60
    time_ticks = []
    for value in range(0, 61, 10):
        x = 75 + 1065 * value / 60
        time_ticks.append(f'<line class="tick" x1="{x:.1f}" y1="350" x2="{x:.1f}" y2="357"/><text class="small" x="{x:.1f}" y="376" text-anchor="middle">{value}</text>')
    fft_ticks = []
    for value in range(0, 7):
        x = 75 + 1065 * value / 6
        fft_ticks.append(f'<line class="tick" x1="{x:.1f}" y1="710" x2="{x:.1f}" y2="717"/><text class="small" x="{x:.1f}" y="736" text-anchor="middle">{value}</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="780" viewBox="0 0 1200 780">
<rect width="1200" height="780" fill="white"/><style>text{{font-family:"Yu Gothic",Meiryo,sans-serif;fill:#163345}}.axis{{stroke:#829aa7;stroke-width:1.2}}.tick{{stroke:#829aa7;stroke-width:1.2}}.grid{{stroke:#e8eff2;stroke-width:1}}.small{{font-size:14px;fill:#607785}}</style>
<text x="600" y="38" text-anchor="middle" font-size="22" font-weight="bold">{html.escape(title)} | 頭部加速度 AP（前後）</text>
<text x="75" y="82" font-size="18" font-weight="bold">時刻波形（全区間）</text><rect x="{shade_x:.1f}" y="105" width="{shade_w:.1f}" height="245" fill="#79c2d0" opacity=".22"/>
<line class="axis" x1="75" y1="350" x2="1140" y2="350"/><line class="axis" x1="75" y1="105" x2="75" y2="350"/>
<polyline points="{wave_pts}" fill="none" stroke="#176b87" stroke-width="1"/>{''.join(time_ticks)}<text class="small" x="607" y="400" text-anchor="middle">Time (s)</text><text class="small" x="1135" y="125" text-anchor="end">色付き: 主解析区間 5–35 s</text>
<text x="75" y="442" font-size="18" font-weight="bold">FFT（主解析区間）</text><line class="axis" x1="75" y1="710" x2="1140" y2="710"/><line class="axis" x1="75" y1="465" x2="75" y2="710"/>
<polyline points="{fft_pts}" fill="none" stroke="#176b87" stroke-width="1.5"/>{''.join(fft_ticks)}<text class="small" x="607" y="770" text-anchor="middle">Frequency (Hz)</text>
<rect x="895" y="480" width="225" height="88" rx="10" fill="#f4f8fa" stroke="#c7d7df"/><text x="1100" y="505" text-anchor="end" font-size="15">RMS: {rms:.3f} mm/s²</text><text x="1100" y="530" text-anchor="end" font-size="15">ピーク: {peak_hz:.2f} Hz</text><text x="1100" y="555" text-anchor="end" font-size="15">全区間RMS: {full_rms:.3f} mm/s²</text>
</svg>'''
    path.write_text(svg, encoding="utf-8")


def comparison_svg(path, records):
    width, height = 1200, 120 + 45 * len(records)
    max_v = max(r[2] for r in records) * 1.12
    rows = []
    for i, (_, stem, rms, _) in enumerate(records):
        y = 92 + i * 45
        bar_w = 700 * rms / max_v
        rows.append(f'<text x="320" y="{y+18}" text-anchor="end" font-size="14">{html.escape(display_name(stem))}</text><rect x="340" y="{y}" width="{bar_w:.1f}" height="26" rx="3" fill="#176b87"/><text x="{350+bar_w:.1f}" y="{y+18}" font-size="14">{rms:.2f}</text>')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:"Yu Gothic",Meiryo,sans-serif;fill:#163345}}</style><text x="600" y="38" text-anchor="middle" font-size="22" font-weight="bold">頭部加速度 AP（前後）RMS 比較</text><text x="600" y="64" text-anchor="middle" font-size="14" fill="#607785">5 ≤ t &lt; 35 s / 平均除去後 RMS (mm/s²)</text>{''.join(rows)}</svg>'''
    path.write_text(svg, encoding="utf-8")


def main():
    OUT.mkdir(exist_ok=True)
    summary = [r for r in read_csv(SOURCE / "rms_summary.csv") if r["axis"] == "AP"]
    summary_fields = list(summary[0])
    write_csv(OUT / "rms_summary_AP.csv", summary_fields, summary)

    records = []
    for row in summary:
        stem = Path(row["file"]).stem
        wave = read_csv(SOURCE / f"{stem}_waveform.csv")
        fft = read_csv(SOURCE / f"{stem}_fft.csv")

        wave_rows = [{
            "time_s": r["time_s"],
            "AP_axhd_mm_s2": r["AP_axhd_mm_s2"],
            "in_5_35s": r["in_5_35s"],
        } for r in wave]
        fft_rows = [{
            "frequency_Hz": r["frequency_Hz"],
            "AP_amplitude_mm_s2": r["AP_amplitude_mm_s2"],
        } for r in fft]
        write_csv(OUT / f"{stem}_waveform_AP.csv", list(wave_rows[0]), wave_rows)
        write_csv(OUT / f"{stem}_fft_AP.csv", list(fft_rows[0]), fft_rows)

        t = [float(r["time_s"]) for r in wave_rows]
        ap = [float(r["AP_axhd_mm_s2"]) for r in wave_rows]
        freq = [float(r["frequency_Hz"]) for r in fft_rows]
        amp = [float(r["AP_amplitude_mm_s2"]) for r in fft_rows]
        rms = float(row["RMS_demeaned_mm_s2"])
        full_rms = float(row["full_0_60s_RMS_demeaned_mm_s2"])
        peak_hz = float(row["FFT_peak_nonDC_Hz"])

        record_svg(OUT / f"{stem}_AP.svg", display_name(stem), t, ap, freq, amp, rms, full_rms, peak_hz)
        records.append((category(stem), stem, rms, peak_hz))

    comparison_svg(OUT / "rms_comparison_AP.svg", records)

    sections = []
    for group in ["静止", "腰部条件", "下肢条件", "移動条件", "その他"]:
        items = [r for r in records if r[0] == group]
        if not items:
            continue
        cards = []
        for _, stem, rms, peak in items:
            cards.append(f'''<article class="card"><h3>{html.escape(display_name(stem))}</h3>
<p class="metrics"><strong>RMS {rms:.3f}</strong> mm/s² <span>ピーク {peak:.2f} Hz</span></p>
<a href="{stem}_AP.svg"><img loading="lazy" src="{stem}_AP.svg" alt="{html.escape(stem)} AP解析図"></a>
<p class="links"><a href="{stem}_waveform_AP.csv">波形CSV</a><a href="{stem}_fft_AP.csv">FFT CSV</a></p></article>''')
        sections.append(f'<section><h2>{group}</h2><div class="grid">{"".join(cards)}</div></section>')

    page = f'''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>頭部IMU AP方向解析</title><style>
:root{{--ink:#163345;--muted:#607785;--line:#dce6eb;--accent:#176b87;--bg:#f4f8fa}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:16px/1.65 "Yu Gothic",Meiryo,sans-serif}}
header,main{{max-width:1180px;margin:auto;padding:28px}}header{{padding-top:48px}}h1{{margin:.1em 0;font-size:clamp(1.8rem,4vw,3rem)}}
.lead{{color:var(--muted);max-width:850px}}.summary{{background:white;border:1px solid var(--line);border-radius:16px;padding:22px;margin:24px 0}}
.summary img{{display:block;width:100%;margin-top:10px}}h2{{margin:48px 0 16px;border-bottom:2px solid var(--accent);padding-bottom:8px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:20px}}.card{{background:white;border:1px solid var(--line);border-radius:16px;padding:18px;box-shadow:0 5px 18px #1730470d}}
.card h3{{margin:0 0 4px}}.card img{{width:100%;display:block;border-radius:9px;border:1px solid var(--line)}}.metrics{{margin:0 0 12px;color:var(--muted)}}
.metrics span{{margin-left:18px}}.links{{display:flex;gap:10px;margin:14px 0 0}}a{{color:var(--accent)}}.links a{{text-decoration:none;border:1px solid var(--line);padding:5px 12px;border-radius:20px}}
.method{{font-size:.94rem;color:var(--muted)}}@media(max-width:560px){{header,main{{padding:18px}}.grid{{grid-template-columns:1fr}}}}
</style></head><body><header><p>HEAD IMU ANALYSIS</p><h1>頭部加速度 AP（前後）方向</h1>
<p class="lead">13記録をAP方向だけに絞り、波形・FFT・RMSを条件別に整理しました。ML（左右）方向の図表と数値は含めていません。</p></header>
<main><section class="summary"><h2>全記録のRMS比較</h2><a href="rms_comparison_AP.svg"><img src="rms_comparison_AP.svg" alt="AP方向RMS比較"></a>
<p><a href="rms_summary_AP.csv">AP方向 RMS一覧CSV</a></p><p class="method">主解析区間: 5 ≤ t &lt; 35 s / 100 Hz / 3000点。RMSは区間平均を除去して算出。FFTは平均除去後にHann窓を適用した片側振幅スペクトルです。</p></section>
{"".join(sections)}</main></body></html>'''
    (OUT / "index.html").write_text(page, encoding="utf-8")
    (OUT / "README.txt").write_text(
        "頭部IMU AP（前後）方向整理版\n\nindex.html を開いてください。\nML（左右）方向は除外しています。\n元データと元フォルダは変更していません。\n",
        encoding="utf-8-sig",
    )
    archive = shutil.make_archive(str(OUT), "zip", root_dir=OUT.parent, base_dir=OUT.name)
    print(f"Created {OUT}")
    print(f"Created {archive}")


if __name__ == "__main__":
    main()
