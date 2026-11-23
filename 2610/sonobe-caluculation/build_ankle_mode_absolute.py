from pathlib import Path

from build_hip_mode_absolute import build


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    build(
        mode_column="xix1",
        strategy_label="足関節",
        device_field="RMSXI1ap",
        out_dir=root / "ankle_strategy_absolute_20261026",
        file_prefix="ankle_mode",
    )
