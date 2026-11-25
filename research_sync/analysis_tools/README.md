# Head IMU gravity compensation (provisional)

## Cycle-resolved visual response analysis

```powershell
.\.venv\Scripts\python.exe -X utf8 research_sync/analysis_tools/analyze_periodic_response.py --output research_sync/results/NEW_OUTPUT
.\.venv\Scripts\python.exe -X utf8 research_sync/analysis_tools/periodic_reference_report.py --output research_sync/results/NEW_OUTPUT
.\.venv\Scripts\python.exe -m unittest discover -s research_sync/analysis_tools -p test_periodic_response.py
```

The first command inventories metadata under `research_sync/recordings` and `2611`, deduplicates trial IDs, and selects complete equal-component-velocity ten-tone trials. It checks logged input frequencies, uses stored clock mapping without optimizing shifts, and estimates native-timestamp harmonic coefficients for four full 20 s cycles. Output CSVs contain per-cycle complex coefficients, aggregated FRFs, coherence, circular phase variation, S/N and input/synchronization audits. Each participant stays separate; figure labels are anonymous. The second command writes a Japanese interpretation and historical rest references and verifies source hashes. Output directories must be new and private (use ignored `research_sync/results`).

The primary available output is provisional AP force/mass COM acceleration. Head Z remains raw sensor counts because lab-coordinate gravity compensation is unvalidated. Sonobe COM position is marked unavailable when participant body inertia/COM height are not confirmed; earlier male/female default curves are not adopted as individual parameters. No original-paper Tv/gamma estimate is reported without independent participant controller parameters and reliable phase. Full assumptions and equations are recorded in the generated RESULTS.md.

## Exploratory head gravity compensation

Run with repository venv dependencies from `research_sync/requirements.txt`:

```powershell
.\.venv\Scripts\python.exe research_sync/analysis_tools/correct_head_gravity.py --output research_sync/results/NEW_LOCAL_OUTPUT PATH_TO_TRIAL_A PATH_TO_TRIAL_B PATH_TO_TRIAL_C
```

Each input directory contains matching `<directory>_metadata.json`, `_quest.csv`, and `_combalance.csv`. Inputs are read only. The output directory must be new; use the ignored `research_sync/results/` tree. Never commit measurements, derived participant data, OEM binaries, or internal inspection dumps.

Plots use anonymous Participant A/B/C labels in argument order. Gray dots show visual displacement, red shows the Sonobe 2019 single-link COM estimate (male coefficients), and blue shows provisional head displacement. Scales are independent. No participant names or local paths are embedded in the script.

External channels 0–2 are acceleration and 3–5 angular velocity, traced in the installed OEM `AdToEng` processing. Scales, zero offsets and range selections come from each recording's device parameters. Only identity axis calibration matrices are currently supported. Initial gyro bias and gravity direction use the two seconds before stimulus, assumed stationary. The initial sensor Z projected onto horizontal defines AP.

Orientation uses gyro integration and accelerometer gravity alignment (30 s time constant). The shaded range compares 10, 30 and 60 s settings; it is NOT a confidence interval. Rotate specific force into the initial frame, subtract gravity, and project horizontally. Displacement uses reflected extension and spectral double integration, removes frequencies below 0.05 Hz, and applies a fourth-order 1 Hz forward-backward low-pass. COM uses Rx=Fy, Ny=-Mx and the same 1 Hz low-pass before selecting 10–90 s; these signs and moment origin remain provisional.

This is not validated head translation. Bias drift, initial motion, unverified gyro calibration, heading drift, translation/tilt ambiguity and boundary assumptions remain. A smaller result is not evidence of accuracy. The setting sensitivity can be large; retain the audit JSON and avoid treating the central curve as a measured position. Original raw values remain available for comparison.
### 前後方向（AP）の基準

今後の床反力・COM解析は、2026年11月15日の mori・takato の計測と同じ被験者・プレートの向きを基準とし、前方を正とする。右踵落とし（両名の cos1／cos2、収録開始約1.4–1.6秒）では4試行とも変換前の `rawMx`、変換後の `My`、前後COP `xp` が負側へ変化した。VRcombalanceの `sensor_Mx_Nm` が同じメーカー座標に対応する記録では、前後COPの符号は `+sensor_Mx_Nm / sensor_Fz_N` を使用する。これはCOPの符号規約であり、モーメント原点の補正や水平床反力項の符号は個別に確認する。配置を変更した場合は変更を明記して符号を再確認する。
