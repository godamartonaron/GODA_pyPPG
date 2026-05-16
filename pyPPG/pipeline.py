"""Programmatic array API for pyPPG.

The original public example is file-oriented and assumes an interactive
workflow.  This module exposes the same preprocessing, fiducial, SQI, and
biomarker steps for callers that already have an in-memory PPG vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

import numpy as np
import pandas as pd
from dotmap import DotMap
from scipy.signal import resample_poly

import pyPPG.preproc as PP
import pyPPG.fiducials as FP
import pyPPG.biomarkers as BM
import pyPPG.ppg_sqi as SQI
from pyPPG import Biomarkers, Fiducials, PPG


HighFsStrategy = Literal["native", "resample", "auto"]
Polarity = Literal["normal", "inverted", 1, -1]


@dataclass
class PPGAnalysisResult:
    """Container returned by :func:`process_ppg_array`."""

    signal: PPG
    fiducials: Fiducials
    fiducials_original_fs: pd.DataFrame
    biomarkers: Biomarkers | None
    ppg_sqi: np.ndarray
    input_fs: float
    analysis_fs: float
    sample_index_scale: float


def _finite_signal(values: np.ndarray) -> np.ndarray:
    x = np.array(values, dtype=np.float64, copy=True).reshape(-1)
    if x.size == 0:
        raise ValueError("PPG signal must not be empty")
    finite = np.isfinite(x)
    if finite.all():
        return x
    if not finite.any():
        raise ValueError("PPG signal contains no finite samples")
    indices = np.arange(x.size)
    x[~finite] = np.interp(indices[~finite], indices[finite], x[finite])
    return x


def _polarity_multiplier(polarity: Polarity) -> int:
    if polarity in ("normal", 1):
        return 1
    if polarity in ("inverted", -1):
        return -1
    raise ValueError("polarity must be 'normal', 'inverted', 1, or -1")


def _resample_if_requested(
    x: np.ndarray,
    fs: float,
    target_fs: float | None,
    strategy: HighFsStrategy,
) -> tuple[np.ndarray, float, float]:
    if strategy not in {"native", "resample", "auto"}:
        raise ValueError("high_fs_strategy must be 'native', 'resample', or 'auto'")
    if fs <= 0:
        raise ValueError("Sampling frequency should be strictly positive")
    if target_fs is None or strategy == "native":
        return x, float(fs), 1.0
    if target_fs <= 0:
        raise ValueError("target_fs should be strictly positive")
    if strategy == "auto" and fs <= target_fs:
        return x, float(fs), 1.0
    if np.isclose(fs, target_fs):
        return x, float(fs), 1.0

    ratio = Fraction(float(target_fs) / float(fs)).limit_denominator(1000)
    y = resample_poly(x, ratio.numerator, ratio.denominator)
    analysis_fs = float(fs) * ratio.numerator / ratio.denominator
    return y, analysis_fs, float(fs) / analysis_fs


def _build_signal(
    values: np.ndarray,
    fs: float,
    *,
    name: str,
    filtering: bool,
    fL: float,
    fH: float,
    order: int,
    sm_wins: dict,
    correction: pd.DataFrame | None,
) -> PPG:
    s = DotMap()
    s.v = values
    s.fs = float(fs)
    s.name = name
    s.start_sig = 0
    s.end_sig = len(values)
    s.filtering = bool(filtering)
    s.fL = fL
    s.fH = fH
    s.order = int(order)
    s.sm_wins = sm_wins
    s.correction = pd.DataFrame() if correction is None else correction.copy()

    prep = PP.Preprocess(fL=fL, fH=fH, order=order, sm_wins=sm_wins)
    s.ppg, s.vpg, s.apg, s.jpg = prep.get_signals(s=s)
    return PPG(s=s, check_ppg_len=True)


def _map_fiducials_to_original_fs(fiducials: pd.DataFrame, scale: float) -> pd.DataFrame:
    mapped = fiducials.copy()
    for column in mapped.columns:
        mapped[column] = np.rint(mapped[column].astype("float64") * scale)
    return mapped.astype("Int64")


def process_ppg_array(
    ppg: np.ndarray,
    fs: float,
    *,
    name: str = "ppg",
    filtering: bool = True,
    fL: float = 0.5000001,
    fH: float = 12,
    order: int = 4,
    sm_wins: dict | None = None,
    correction: pd.DataFrame | None = None,
    high_fs_strategy: HighFsStrategy = "native",
    target_fs: float | None = None,
    polarity: Polarity = "normal",
    get_biomarkers: bool = True,
) -> PPGAnalysisResult:
    """Run pyPPG on an in-memory PPG vector.

    High-frequency PPG is supported in two ways.  With the default
    ``high_fs_strategy="native"``, filtering uses numerically stable SOS
    filters at the original sampling frequency.  With ``"resample"`` or
    ``"auto"``, ``target_fs`` can be used to analyze a downsampled signal and
    ``fiducials_original_fs`` maps detected sample indices back to the input
    sampling frequency.
    """

    sm_wins = {"ppg": 50, "vpg": 10, "apg": 10, "jpg": 10} if sm_wins is None else sm_wins
    raw = _finite_signal(ppg) * _polarity_multiplier(polarity)
    analysis_signal, analysis_fs, sample_index_scale = _resample_if_requested(
        raw,
        float(fs),
        target_fs,
        high_fs_strategy,
    )

    signal = _build_signal(
        analysis_signal,
        analysis_fs,
        name=name,
        filtering=filtering,
        fL=fL,
        fH=fH,
        order=order,
        sm_wins=sm_wins,
        correction=correction,
    )

    fpex = FP.FpCollection(s=signal)
    fiducial_frame = fpex.get_fiducials(s=signal)
    fiducials = Fiducials(fp=fiducial_frame)
    fiducials_original_fs = _map_fiducials_to_original_fs(fiducial_frame, sample_index_scale)

    sp = fiducial_frame["sp"].dropna().astype(int).to_numpy()
    ppg_sqi = SQI.get_ppgSQI(ppg=signal.ppg, fs=signal.fs, annotation=sp)

    biomarkers = None
    if get_biomarkers and len(fiducial_frame) > 0:
        bmex = BM.BmCollection(s=signal, fp=fiducials)
        bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()
        biomarkers = Biomarkers(bm_defs=bm_defs, bm_vals=bm_vals, bm_stats=bm_stats)

    return PPGAnalysisResult(
        signal=signal,
        fiducials=fiducials,
        fiducials_original_fs=fiducials_original_fs,
        biomarkers=biomarkers,
        ppg_sqi=ppg_sqi,
        input_fs=float(fs),
        analysis_fs=analysis_fs,
        sample_index_scale=sample_index_scale,
    )
