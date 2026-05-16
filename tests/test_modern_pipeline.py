import numpy as np
import pandas as pd

import pyPPG


def synthetic_ppg(fs=1000.0, duration=35.0, heart_rate_hz=1.2):
    t = np.arange(0, duration, 1 / fs)
    primary = np.clip(np.sin(2 * np.pi * heart_rate_hz * t), 0, None) ** 2
    reflection = 0.2 * np.roll(primary, int(0.20 * fs))
    noise = 0.01 * np.random.default_rng(42).normal(size=t.size)
    return primary + reflection + noise


def test_native_high_sampling_rate_pipeline_is_finite():
    fs = 1000.0
    result = pyPPG.process_ppg_array(
        synthetic_ppg(fs=fs),
        fs,
        fH=12,
        high_fs_strategy="native",
        get_biomarkers=True,
    )

    fiducials = result.fiducials.get_fp()
    assert result.analysis_fs == fs
    assert result.sample_index_scale == 1.0
    assert np.isfinite(result.signal.ppg).all()
    assert np.isfinite(result.signal.vpg).all()
    assert len(fiducials) >= 20
    assert fiducials["sp"].notna().all()
    assert result.ppg_sqi.size == len(fiducials) - 1
    assert result.biomarkers is not None
    assert "IPA" in result.biomarkers.bm_stats["sig_ratios"].columns


def test_resampled_high_sampling_rate_maps_fiducials_to_original_indices():
    fs = 1000.0
    result = pyPPG.process_ppg_array(
        synthetic_ppg(fs=fs),
        fs,
        fH=12,
        high_fs_strategy="resample",
        target_fs=125,
        get_biomarkers=False,
    )

    fiducials = result.fiducials.get_fp()
    assert result.analysis_fs == 125.0
    assert result.sample_index_scale == 8.0
    assert len(fiducials) >= 20
    np.testing.assert_allclose(
        result.fiducials_original_fs["sp"].astype("float64").to_numpy(),
        np.rint(fiducials["sp"].astype("float64").to_numpy() * 8.0),
    )


def test_polarity_option_handles_inverted_ppg():
    fs = 1000.0
    result = pyPPG.process_ppg_array(
        -synthetic_ppg(fs=fs),
        fs,
        fH=12,
        polarity="inverted",
        get_biomarkers=False,
    )

    fiducials = result.fiducials.get_fp()
    assert len(fiducials) >= 20
    assert fiducials["sp"].notna().all()


def test_pandas_copy_on_write_keeps_fiducial_values():
    fs = 125.0
    result = pyPPG.process_ppg_array(
        synthetic_ppg(fs=fs),
        fs,
        high_fs_strategy="native",
        get_biomarkers=False,
    )

    fiducials = result.fiducials.get_fp()
    assert isinstance(fiducials, pd.DataFrame)
    assert fiducials[["on", "sp", "off"]].notna().all().all()


def test_ppg_sqi_accepts_nullable_annotation_values():
    from pyPPG.ppg_sqi import get_ppgSQI

    fs = 125.0
    signal = synthetic_ppg(fs=fs, duration=8)
    annotation = pd.Series([25, pd.NA, 125, 225, np.nan, 325], dtype="Float64")

    sqi = get_ppgSQI(signal, fs, annotation)

    assert sqi.shape == (3,)
    assert np.isfinite(sqi).any() or np.isnan(sqi).all()


def test_process_ppg_array_does_not_mutate_input_signal():
    fs = 1000.0
    signal = synthetic_ppg(fs=fs)
    signal[10] = np.nan
    original = signal.copy()

    pyPPG.process_ppg_array(signal, fs, get_biomarkers=False)

    np.testing.assert_equal(signal, original)
