import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = "C:/Users/seguchi/gitseguchi/All/2610/outputs/fft_1027_20261101";
const data = JSON.parse(await fs.readFile(`${outputDir}/fft_data.json`, "utf8"));
const results = data.results;
const workbook = Workbook.create();
const sheet = workbook.worksheets.add("FFT結果");
sheet.showGridLines = false;
sheet.tabColor = "#1F4E78";

const font = "Arial";
sheet.getRange("A1:P40").format.font = { name: font, size: 10, color: "#1F2937" };
sheet.getRange("A2").values = [["頭部加速度 z方向 FFT（5〜30秒、0〜1 Hz）"]];
sheet.getRange("A2").format.font = { name: font, size: 15, bold: true, color: "#1F2937" };
sheet.getRange("A3:P3").format.borders = { bottom: { style: "thin", color: "#94A3B8" } };

const frequencies = results[0].frequencies_hz;
const headers = ["周波数 (Hz)", ...results.map((item) => item.condition)];
const rows = frequencies.map((frequency, index) => [
  frequency,
  ...results.map((item) => item.amplitudes[index]),
]);
sheet.getRange("A4:E4").values = [headers];
sheet.getRange("A5").write(rows);
sheet.getRange("A4:E4").format = {
  fill: "#1F4E78",
  font: { name: font, size: 10, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
  verticalAlignment: "center",
  borders: { insideVertical: { style: "thin", color: "#FFFFFF" } },
};
sheet.getRange(`A5:A${4 + rows.length}`).format.numberFormat = "0.00";
sheet.getRange(`B5:E${4 + rows.length}`).format.numberFormat = "0.000";
sheet.getRange(`A5:E${4 + rows.length}`).format.verticalAlignment = "center";
sheet.getRange(`A5:E${4 + rows.length}`).format.borders = {
  insideHorizontal: { style: "thin", color: "#E2E8F0" },
};

const chart = sheet.charts.add("line", sheet.getRange(`A4:E${4 + rows.length}`));
chart.title = "頭部加速度 z方向の振幅スペクトル";
chart.titleTextStyle.fontSize = 12;
chart.titleTextStyle.typeface = font;
chart.legend = { position: "top", textStyle: { typeface: font, fontSize: 10 } };
chart.xAxis = {
  axisType: "textAxis",
  textStyle: { typeface: font, fontSize: 9 },
  numberFormatCode: "0.00",
  numberFormatSourceLinked: false,
};
chart.yAxis = {
  textStyle: { typeface: font, fontSize: 9 },
  numberFormatCode: "0.0",
  numberFormatSourceLinked: false,
};
chart.xAxis.title.text = "周波数 (Hz)";
chart.yAxis.title.text = "片側振幅（元データ単位）";
chart.setPosition("G4", "P23");
const colors = ["#2563EB", "#DC2626", "#16A34A", "#9333EA"];
chart.series.items.forEach((series, index) => {
  series.line = { fill: colors[index], style: "solid", width: 2 };
});

sheet.getRange("A32:E32").values = [["処理条件", "値", "単位", "備考", ""]];
sheet.getRange("A32:E32").format = {
  fill: "#DCE6F1",
  font: { name: font, size: 10, bold: true, color: "#1F2937" },
  horizontalAlignment: "center",
};
sheet.getRange("A33:D37").values = [
  ["解析区間", "5.00〜30.00", "秒", "先頭時刻を0秒とした25秒間"],
  ["サンプリング周波数", results[0].sample_rate_hz, "Hz", "時刻差の中央値から算出"],
  ["サンプル数", results[0].window_samples, "点", "各条件共通"],
  ["周波数分解能", results[0].frequency_resolution_hz, "Hz", "25秒窓のため0.04 Hz"],
  ["FFT処理", "平均値除去・片側振幅", "", "窓関数なし"],
];
sheet.getRange("B34").format.numberFormat = "0.00";
sheet.getRange("B35").format.numberFormat = "0";
sheet.getRange("B36").format.numberFormat = "0.00";
sheet.getRange("A39:G39").values = [["条件", "総サンプル数", "解析サンプル数", "解析終端 (秒)", "元ファイル", "", ""]];
sheet.getRange("A39:G39").format = {
  fill: "#1F4E78",
  font: { name: font, size: 10, bold: true, color: "#FFFFFF" },
  horizontalAlignment: "center",
};
sheet.getRange("A40").write(results.map((item) => [
  item.condition,
  item.total_samples,
  item.window_samples,
  item.window_end_seconds,
  item.source_file,
]));
sheet.getRange(`D40:D${39 + results.length}`).format.numberFormat = "0.00";

sheet.getRange("A:A").format.columnWidth = 16;
sheet.getRange("B:E").format.columnWidth = 17;
sheet.getRange("F:F").format.columnWidth = 2;
sheet.getRange("G:P").format.columnWidth = 11;
sheet.getRange("D33:D37").format.columnWidth = 32;
sheet.getRange("E40:E43").format.columnWidth = 58;
sheet.getRange("A2:P43").format.verticalAlignment = "center";

workbook.recalculate();
const inspected = await workbook.inspect({
  kind: "table",
  range: "FFT結果!A2:E43",
  include: "values,formulas",
  tableMaxRows: 45,
  tableMaxCols: 5,
  maxChars: 12000,
});
console.log(inspected.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);
const preview = await workbook.render({ sheetName: "FFT結果", range: "A1:P43", scale: 1.5, format: "png" });
await fs.writeFile(`${outputDir}/頭部加速度z_FFT_グラフ.png`, new Uint8Array(await preview.arrayBuffer()));
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(`${outputDir}/頭部加速度z_FFT_5-30秒_0-1Hz.xlsx`);
