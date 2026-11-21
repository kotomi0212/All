from pathlib import Path
import csv
import html
import io
import math
import shutil
import zipfile


SOURCE_ZIP = Path(r"C:\Users\seguchi\Downloads\1022結果-20261026T091514Z-1-001.zip")
OUT = Path(__file__).resolve().parent / "hip_strategy_absolute_20261026"


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def fft_amplitude_0_6(values, fs=100.0):
    n = len(values)
    mean = sum(values) / n
    centered = [v - mean for v in values]
    window = [0.5 - 0.5 * math.cos(2 * math.pi * i / (n - 1)) for i in range(n)]
    wsum = sum(window)
    result = []
    max_bin = int(6 * n / fs)
    for k in range(max_bin + 1):
        re = im = 0.0
        factor = -2 * math.pi * k / n
        for i, (x, w) in enumerate(zip(centered, window)):
            angle = factor * i
            xw = x * w
            re += xw * math.cos(angle)
            im += xw * math.sin(angle)
        amp = math.hypot(re, im) / wsum
        if k != 0:
            amp *= 2
        result.append((k * fs / n, amp))
    return result


def polyline(values, x0, y0, width, height, ymin=None, ymax=None):
    if ymin is None: ymin = min(values)
    if ymax is None: ymax = max(values)
    span = ymax - ymin or 1
    step = max(1, len(values) // 1800)
    return " ".join(
        f"{x0 + width*i/max(1,len(values)-1):.1f},{y0 + height*(ymax-values[i])/span:.1f}"
        for i in range(0, len(values), step)
    )


def nice_name(stem):
    return stem.removeprefix("261022_").replace("_seguchi", " / ")


def group_name(stem):
    s = stem.lower()
    if "statics" in s: return "静止"
    if "hip" in s: return "股関節中心回転"
    if "leg" in s: return "足関節中心回転"
    if "trans" in s: return "並進"
    return "その他"


def make_svg(path, title, time, signal, spectrum, rms, device_rms, ratio,
             strategy_label="股関節", mode_column="xix2", device_field="RMSXI2ap"):
    freqs = [r[0] for r in spectrum]
    amps = [r[1] for r in spectrum]
    wave = polyline(signal, 75, 105, 1065, 245)
    fft = polyline(amps, 75, 465, 1065, 245, 0, max(amps) * 1.08)
    time_ticks = []
    for v in range(0, 61, 10):
        x = 75 + 1065*v/60
        time_ticks.append(f'<line class="tick" x1="{x}" y1="350" x2="{x}" y2="357"/><text class="small" x="{x}" y="376" text-anchor="middle">{v}</text>')
    fft_ticks = []
    for v in range(7):
        x = 75 + 1065*v/6
        fft_ticks.append(f'<line class="tick" x1="{x}" y1="710" x2="{x}" y2="717"/><text class="small" x="{x}" y="736" text-anchor="middle">{v}</text>')
    shade_x, shade_w = 75 + 1065*5/60, 1065*30/60
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="790" viewBox="0 0 1200 790"><rect width="100%" height="100%" fill="white"/>
<style>text{{font-family:"Yu Gothic",Meiryo,sans-serif;fill:#163345}}.axis,.tick{{stroke:#829aa7;stroke-width:1.2}}.small{{font-size:14px;fill:#607785}}</style>
<text x="600" y="38" text-anchor="middle" font-size="22" font-weight="bold">{html.escape(title)} | {strategy_label}戦略モード絶対量（AP）</text>
<text x="75" y="82" font-size="18" font-weight="bold">時刻波形 {mode_column}（全区間）</text><rect x="{shade_x}" y="105" width="{shade_w}" height="245" fill="#d79a35" opacity=".2"/>
<line class="axis" x1="75" y1="350" x2="1140" y2="350"/><line class="axis" x1="75" y1="105" x2="75" y2="350"/><polyline points="{wave}" fill="none" stroke="#b86420" stroke-width="1"/>{''.join(time_ticks)}
<text class="small" x="607" y="400" text-anchor="middle">Time (s)</text><text class="small" x="1135" y="125" text-anchor="end">色付き: 主解析区間 5–35 s</text>
<text x="75" y="442" font-size="18" font-weight="bold">FFT（5–35 s）</text><line class="axis" x1="75" y1="710" x2="1140" y2="710"/><line class="axis" x1="75" y1="465" x2="75" y2="710"/><polyline points="{fft}" fill="none" stroke="#b86420" stroke-width="1.5"/>{''.join(fft_ticks)}
<text class="small" x="607" y="770" text-anchor="middle">Frequency (Hz)</text><rect x="870" y="478" width="250" height="110" rx="10" fill="#fff8ec" stroke="#e0c28f"/>
<text x="1100" y="505" text-anchor="end" font-size="15">解析区間RMS: {rms:.3f} mm/s²</text><text x="1100" y="532" text-anchor="end" font-size="15">装置{device_field}: {device_rms:.3f} mm/s²</text><text x="1100" y="559" text-anchor="end" font-size="15">RXIap (XI2/XI1): {ratio:.3f}</text>
</svg>'''
    path.write_text(svg, encoding="utf-8")


def comparison_svg(path, records, strategy_label="股関節", mode_column="xix2"):
    h = 120 + 45*len(records)
    max_v = max(r[2] for r in records)*1.12
    rows=[]
    for i, (_, stem, rms, _, _) in enumerate(records):
        y=92+i*45; w=700*rms/max_v
        rows.append(f'<text x="320" y="{y+18}" text-anchor="end" font-size="14">{html.escape(nice_name(stem))}</text><rect x="340" y="{y}" width="{w:.1f}" height="26" rx="3" fill="#b86420"/><text x="{350+w:.1f}" y="{y+18}" font-size="14">{rms:.2f}</text>')
    path.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{h}" viewBox="0 0 1200 {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:"Yu Gothic",Meiryo,sans-serif;fill:#163345}}</style><text x="600" y="38" text-anchor="middle" font-size="22" font-weight="bold">{strategy_label}戦略モード絶対量 RMS 比較</text><text x="600" y="64" text-anchor="middle" font-size="14">{mode_column} / AP / 5 ≤ t &lt; 35 s / mm/s²</text>{''.join(rows)}</svg>''', encoding="utf-8")


def build(mode_column="xix2", strategy_label="股関節", device_field="RMSXI2ap",
          out_dir=OUT, file_prefix="hip_mode"):
    out_dir.mkdir(exist_ok=True)
    records=[]; summary=[]
    with zipfile.ZipFile(SOURCE_ZIP) as z:
        names=sorted(n for n in z.namelist() if n.lower().endswith('.csv'))
        for member in names:
            raw=z.read(member)
            try: text=raw.decode('utf-8-sig')
            except UnicodeDecodeError: text=raw.decode('cp932')
            rows=list(csv.reader(io.StringIO(text)))
            k=max(i for i,r in enumerate(rows) if 'xix2' in r and 'axhd' in r)
            header=rows[k]
            ti=0; xi=header.index(mode_column)
            data=[r for r in rows[k+1:] if len(r)>xi and r[0].strip()]
            time=[float(r[ti]) for r in data]; x=[float(r[xi]) for r in data]
            selected=[v for t,v in zip(time,x) if 5-1e-9 <= t < 35-1e-9]
            mean=sum(selected)/len(selected)
            rms=math.sqrt(sum((v-mean)**2 for v in selected)/len(selected))
            spectrum=fft_amplitude_0_6(selected)
            peak=max(spectrum[1:], key=lambda p:p[1])
            eval_header=next(r for r in rows[:k] if device_field in r)
            eval_row=rows[rows.index(eval_header)+2]
            device_rms=float(eval_row[eval_header.index(device_field)])
            ratio=float(eval_row[eval_header.index('RXIap')])
            stem=Path(member).stem
            write_csv(out_dir/f'{stem}_{file_prefix}_waveform.csv',['time_s',f'{mode_column}_mm_s2','in_5_35s'],zip(time,x,[int(5<=t<35) for t in time]))
            write_csv(out_dir/f'{stem}_{file_prefix}_fft.csv',['frequency_Hz',f'{mode_column}_amplitude_mm_s2'],spectrum)
            make_svg(out_dir/f'{stem}_{file_prefix}.svg',nice_name(stem),time,x,spectrum,rms,device_rms,ratio,strategy_label,mode_column,device_field)
            records.append((group_name(stem),stem,rms,device_rms,peak))
            summary.append([stem,rms,device_rms,ratio,peak[0],peak[1]])
    write_csv(out_dir/f'{file_prefix}_summary.csv',['record',f'RMS_{mode_column}_5_35s_demeaned_mm_s2',f'device_{device_field}_mm_s2','RXIap_ratio_XI2_over_XI1','FFT_peak_0_6Hz','FFT_peak_amplitude_mm_s2'],summary)
    comparison_svg(out_dir/f'{file_prefix}_rms_comparison.svg',records,strategy_label,mode_column)
    sections=[]
    for group in ['静止','股関節中心回転','足関節中心回転','並進','その他']:
        items=[r for r in records if r[0]==group]
        if not items: continue
        cards=[]
        for _,stem,rms,_,peak in items:
            cards.append(f'''<article><h3>{html.escape(nice_name(stem))}</h3><p><strong>RMS {rms:.3f}</strong> mm/s² / ピーク {peak[0]:.2f} Hz</p><a href="{stem}_{file_prefix}.svg"><img src="{stem}_{file_prefix}.svg"></a><p><a href="{stem}_{file_prefix}_waveform.csv">波形CSV</a> <a href="{stem}_{file_prefix}_fft.csv">FFT CSV</a></p></article>''')
        sections.append(f'<section><h2>{group}</h2><div class="grid">{"".join(cards)}</div></section>')
    page=f'''<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{strategy_label}戦略モード絶対量</title><style>body{{margin:0;background:#f7f5f1;color:#163345;font:16px/1.65 "Yu Gothic",Meiryo,sans-serif}}header,main{{max-width:1180px;margin:auto;padding:28px}}h1{{font-size:clamp(1.8rem,4vw,3rem)}}.note,article{{background:white;border:1px solid #e6ddd0;border-radius:16px;padding:20px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:20px}}img{{width:100%;display:block}}h2{{margin-top:48px;border-bottom:2px solid #b86420;padding-bottom:8px}}a{{color:#9b4d17;margin-right:12px}}@media(max-width:560px){{header,main{{padding:18px}}.grid{{grid-template-columns:1fr}}}}</style></head><body><header><h1>{strategy_label}戦略の絶対的な揺れ</h1><p>寄与率ではなく、割合化前の{strategy_label}戦略モード <code>{mode_column}</code>（AP方向、mm/s²）を解析しています。これは関節角度ではなく、身体運動を戦略基底へ分解したモード加速度相当量です。</p></header><main><section class="note"><h2>全記録のRMS比較</h2><img src="{file_prefix}_rms_comparison.svg"><p><a href="{file_prefix}_summary.csv">集計CSV</a></p><p>FFT: 5–35秒、平均除去、Hann窓、片側振幅、0–6 Hz。横軸は1 Hz刻みです。</p></section>{''.join(sections)}</main></body></html>'''
    (out_dir/'index.html').write_text(page,encoding='utf-8')
    (out_dir/'README.txt').write_text(f'{strategy_label}戦略モード絶対量の解析版。{mode_column}（AP, mm/s²）を割合化せずにFFT・RMS解析。\n',encoding='utf-8-sig')
    shutil.make_archive(str(out_dir),'zip',root_dir=out_dir.parent,base_dir=out_dir.name)
    print(f'Created {out_dir} ({len(records)} records)')


def main():
    build()


if __name__=='__main__': main()
