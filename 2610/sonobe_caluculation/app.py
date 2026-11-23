"""Local desktop viewer and command-line entry point."""
import argparse
from pathlib import Path
import sys

from balance_analysis import read_recording, analyze, make_figure, export_results


def launch(initial=None):
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    root = tk.Tk()
    root.title('立位バランス | 足関節・股関節モード解析')
    root.geometry('1240x960')
    root.minsize(960, 760)
    style = ttk.Style(root)
    style.theme_use('clam')
    style.configure('.', font=('Yu Gothic UI', 10))
    state = {'recording': None, 'result': None, 'canvas': None}
    header = ttk.Frame(root, padding=(18, 12))
    header.pack(fill='x')
    ttk.Label(header, text='足関節・股関節モード解析', font=('Yu Gothic UI', 19, 'bold')).pack(anchor='w')
    ttk.Label(header, text='Sonobe & Mitsutake (2025)  ·  MedicalForcePlate 処理済みCSV').pack(anchor='w')
    controls = ttk.Frame(root, padding=(18, 4))
    controls.pack(fill='x')
    start, end = tk.StringVar(value='5'), tk.StringVar(value='35')
    status = tk.StringVar(value='CSVを選択してください。解析はこのPC内で実行します。')
    note = tk.StringVar(value='割合は各モードRMS / 両モードRMS合計。時間やトルクの割合ではありません。')
    ttk.Label(root, textvariable=status, padding=(18, 6), wraplength=1150).pack(fill='x')
    plot_frame = ttk.Frame(root)
    plot_frame.pack(fill='both', expand=True)
    toolbar_frame = ttk.Frame(root)
    toolbar_frame.pack(fill='x')
    ttk.Label(root, textvariable=note, padding=(18, 10), wraplength=1150,
              foreground='#844710').pack(fill='x')
    ttk.Label(root, text='平均を除去してRMSを計算。追加フィルタなし。MLの股関節モードは腰部モデルに対応。',
              padding=(18, 4)).pack(fill='x')

    def calculate():
        if state['recording'] is None:
            return
        try:
            result = analyze(state['recording'], float(start.get()), float(end.get()))
        except ValueError as exc:
            state['result'] = None
            for widget in plot_frame.winfo_children():
                widget.destroy()
            for widget in toolbar_frame.winfo_children():
                widget.destroy()
            note.set('解析できません。区間を修正して再計算してください。')
            messagebox.showerror('解析条件を確認してください', str(exc))
            return
        for widget in plot_frame.winfo_children():
            widget.destroy()
        for widget in toolbar_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(make_figure(state['recording'], result), master=plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)
        NavigationToolbar2Tk(canvas, toolbar_frame)
        state.update(result=result, canvas=canvas)
        r = state['recording']
        status.set(f'{r.path.name}  |  {r.fs:g} Hz  |  {len(r.time):,}点  |  '
                   f"解析: {result['start_s']:g} ≤ t < {result['end_s']:g} 秒  ({result['samples']:,}点)")
        note.set('\n'.join(result['warnings']) or 'モード和と頭部加速度の整合性検査：閾値内（頭部RMSの1%）。')

    def open_file(path=None):
        path = path or filedialog.askopenfilename(filetypes=[('CSV', '*.csv'), ('All files', '*.*')])
        if not path:
            return
        try:
            recording = read_recording(path)
        except (ValueError, OSError) as exc:
            messagebox.showerror('読み込みできません', str(exc))
            return
        state['recording'] = recording
        state['result'] = None
        # Do not silently change the standard 5–35 s interval for short files.
        start.set('5')
        end.set('35')
        if recording.time[0] > 5 or recording.time[-1]+1/recording.fs < 35-1e-7:
            start.set(f'{recording.time[0]:g}')
            end.set(f'{recording.time[-1]+1/recording.fs:g}')
            messagebox.showinfo('解析区間', '5〜35秒が記録範囲外なので、全区間を設定しました。画面の区間を確認してください。')
        calculate()

    def full_range():
        r = state['recording']
        if r is not None:
            start.set(f'{r.time[0]:g}')
            end.set(f'{r.time[-1]+1/r.fs:g}')
            calculate()

    def save():
        if state['result'] is None:
            messagebox.showinfo('結果保存', '先にCSVを読み込み、解析してください。')
            return
        destination = filedialog.askdirectory(title='結果フォルダーの保存先を選択')
        if destination:
            try:
                out = export_results(state['recording'], state['result'], destination)
                messagebox.showinfo('保存完了', f'表示中の解析結果を保存しました。\n{out}\n\nPNG・指標CSV・波形CSV・JSON・注意事項')
            except OSError as exc:
                messagebox.showerror('保存できません', str(exc))

    ttk.Button(controls, text='CSVを開く', command=open_file).pack(side='left', padx=(0, 20))
    ttk.Label(controls, text='開始 (秒)').pack(side='left')
    ttk.Entry(controls, textvariable=start, width=7).pack(side='left', padx=5)
    ttk.Label(controls, text='終了 (秒・含まない)').pack(side='left')
    ttk.Entry(controls, textvariable=end, width=7).pack(side='left', padx=5)
    ttk.Button(controls, text='再計算', command=calculate).pack(side='left', padx=8)
    ttk.Button(controls, text='全区間', command=full_range).pack(side='left')
    ttk.Button(controls, text='表示中の結果を保存', command=save).pack(side='right')
    if initial:
        root.after(100, lambda: open_file(initial))
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description='Sonobe 2025: 処理済みモードCSVの割合解析')
    parser.add_argument('csv', nargs='?', type=Path)
    parser.add_argument('--start', type=float, default=5)
    parser.add_argument('--end', type=float, default=35)
    parser.add_argument('--output', type=Path, help='GUIなしで結果を書き出す保存先')
    args = parser.parse_args()
    if args.output:
        if args.csv is None:
            parser.error('--output にはCSVパスが必要です。')
        try:
            recording = read_recording(args.csv)
            result = analyze(recording, args.start, args.end)
            target = export_results(recording, result, args.output)
        except (ValueError, OSError) as exc:
            parser.exit(1, f'Error: {exc}\n')
        print(target)
        for axis, m in result['axes'].items():
            print(axis, m)
        for warning in result['warnings']:
            print('注意:', warning)
    else:
        launch(args.csv)


if __name__ == '__main__':
    main()
