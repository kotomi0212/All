# Head IMU gravity compensation (provisional)

Run with repository venv dependencies from `research_sync/requirements.txt`:

```powershell
.\.venv\Scripts\python.exe research_sync/analysis_tools/correct_head_gravity.py --output research_sync/results/NEW_LOCAL_OUTPUT PATH_TO_TRIAL_A PATH_TO_TRIAL_B PATH_TO_TRIAL_C
```

Each input directory contains matching `<directory>_metadata.json`, `_quest.csv`, and `_combalance.csv`. Inputs are read only. The output directory must be new; use the ignored `research_sync/results/` tree. Never commit measurements, derived participant data, OEM binaries, or internal inspection dumps.

Plots use anonymous Participant A/B/C labels in argument order. Gray dots show visual displacement, red shows the Sonobe 2019 single-link COM estimate (male coefficients), and blue shows provisional head displacement. Scales are independent. No participant names or local paths are embedded in the script.

External channels 0–2 are acceleration and 3–5 angular velocity, traced in the installed OEM `AdToEng` processing. Scales, zero offsets and range selections come from each recording's device parameters. Only identity axis calibration matrices are currently supported. Initial gyro bias and gravity direction use the two seconds before stimulus, assumed stationary. The initial sensor Z projected onto horizontal defines AP.

Orientation uses gyro integration and accelerometer gravity alignment (30 s time constant). The shaded range compares 10, 30 and 60 s settings; it is NOT a confidence interval. Rotate specific force into the initial frame, subtract gravity, and project horizontally. Displacement uses reflected extension and spectral double integration, removes frequencies below 0.05 Hz, and applies a fourth-order 1 Hz forward-backward low-pass. COM uses Rx=Fy, Ny=-Mx and the same 1 Hz low-pass before selecting 10–90 s; these signs and moment origin remain provisional.

This is not validated head translation. Bias drift, initial motion, unverified gyro calibration, heading drift, translation/tilt ambiguity and boundary assumptions remain. A smaller result is not evidence of accuracy. The setting sensitivity can be large; retain the audit JSON and avoid treating the central curve as a measured position. Original raw values remain available for comparison.
