from pathlib import Path

import pandas as pd


def test_stage_mapping_file_exists() -> None:
    mapping = Path("code/coarse_graining/anatomical7_mapping_from_dmt.csv")
    assert mapping.exists()


def test_stage_labels_are_standardized() -> None:
    cfg_path = Path("code/run_all/somnodynamics_config.toml")
    assert cfg_path.exists()
    text = cfg_path.read_text(encoding="utf-8")
    assert 'W_ = "Wake"' in text
    assert 'R_ = "REM"' in text


def test_macro_regions_are_expected() -> None:
    df = pd.read_csv("code/coarse_graining/anatomical7_mapping_from_dmt.csv")
    expected = {
        "Frontal",
        "Sensorimotor",
        "Parietal",
        "Temporal",
        "Occipital",
        "Cingulate",
        "Limbic",
    }
    assert set(df["macro_region"]).issubset(expected)
