# 視覚刺激 → COM加速度の自動同定

MATLAB R2025b + System Identification Toolbox向け。既存の
`2611/261121/fit_visual_acceleration_first_order.py` のCSV列・時計対応・
AP方向（90度反時計回りの立位でFy）を参照した。

MATLABでこのフォルダをパスに追加し、実行する。

```matlab
result = analyzeVisualTrial; % 試行フォルダ選択
% または
result = analyzeVisualTrial("C:\Users\seguchi\gitseguchi\All\2611\261121\20261121_otsubo_90s_multicosine_001");
```

入力は実際に適用した `ap_acceleration_m_s2`、出力は `sensor_Fy_N / body_mass_kg`。
PC対応時刻で同期し、100 Hzの等間隔データへ線形補間する。
既定区間10–90秒を前半60%の推定区間と後半40%の検証区間に分ける。
`procest` と `idproc('P1')` により `G(s)=-gamma/(Tv*s+1)` を推定する。
ゲインの符号を負、時定数を正に制約する。遅延は推定しない。
既存Pythonの周波数領域・対数ゲインと位相の目的関数とは異なるため、
推定値が一致することは要求しない。帯域制限は行わず時間波形を同定する。

新規の結果フォルダにMATモデル、指標CSV、同期波形CSV、検証図、残差図を保存する。
`result.Model` と `result.Data` はSystem Identificationアプリへのインポートにも利用できる。

```matlab
% 別刺激の解析区間、実測体重、立位方向を明示する場合
result = analyzeVisualTrial(folder,Window=[10 90],MassKg=62.0, ...
    ForceColumn="sensor_Fy_N",ForceSign=1);
```

単一刺激では試行の定常区間に合わせてWindowを変更する。
単一周波数だけでTvとgammaの再現性を確定しない。
サンプリングを粗くする場合は別途アンチエイリアス処理が必要。
欠損、時計の逆転、100 ms以上の記録断、無入力は解析を停止する。

品質判定は暫定。検証適合率50%超、ゼロ応答より小さいRMSE、時定数の
境界張り付きなしの場合のみPROVISIONAL_FIT、それ以外はREVIEW_REQUIRED。
これはS/N比や生理学的妥当性の証明ではない。残差、入力周波数の応答、
静止ノイズ、別試行での再現性を合わせて確認する。
ログの体重は未検証の静止Fz推定値の場合がある。時計対応も未校正のため、
位相とTvの解釈には計測遅延の確認が必要。
