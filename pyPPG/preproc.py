import numpy as np
from dotmap import DotMap
from scipy import signal


def _as_finite_1d(values):
    x = np.asarray(values, dtype=np.float64).reshape(-1)
    if x.size == 0:
        return x
    finite = np.isfinite(x)
    if finite.all():
        return x
    if not finite.any():
        return np.zeros_like(x)
    indices = np.arange(x.size)
    x[~finite] = np.interp(indices[~finite], indices[finite], x[finite])
    return x


def _safe_cutoff(value, fs, *, low=False):
    nyquist = float(fs) / 2.0
    if nyquist <= 0:
        raise ValueError("Sampling frequency should be strictly positive")
    cutoff = float(value)
    if low:
        return max(0.0, min(cutoff, nyquist * 0.95))
    return max(1e-6, min(cutoff, nyquist * 0.95))


def _zero_phase_sos_filter(x, fs, fL, fH, order):
    x = _as_finite_1d(x)
    if x.size < 4:
        return x

    high = _safe_cutoff(fH, fs)
    if fL == 0:
        sos = signal.cheby2(order, 20, high, btype="lowpass", fs=fs, output="sos")
    else:
        low = _safe_cutoff(fL, fs, low=True)
        if high <= low:
            return x - np.nanmedian(x)
        sos = signal.cheby2(order, 20, [low, high], btype="bandpass", fs=fs, output="sos")

    return signal.sosfiltfilt(sos, x)


def _window_samples(fs, milliseconds):
    win = int(round(float(fs) * float(milliseconds) / 1000.0))
    return max(1, win)


def _moving_average_zero_phase(x, win):
    x = _as_finite_1d(x)
    if win <= 1 or x.size < 3:
        return x
    kernel = np.ones(win, dtype=np.float64) / float(win)
    padlen = 3 * (max(len(kernel), 1) - 1)
    if x.size <= padlen:
        return np.convolve(x, kernel, mode="same")
    return signal.filtfilt(kernel, [1.0], x)


class Preprocess:

    ###########################################################################
    ######################## Initialization of Biomarkers #####################
    ###########################################################################
    def __init__(self,fL=0.5000001, fH=12, order=4, sm_wins={'ppg':50,'vpg':10,'apg':10,'jpg':10}):
        """
        The purpose of the Preprocess class is to filter and calculate the PPG, PPG', PPG", and PPG'" signals.

        :param fL: Lower cutoff frequency (Hz)
        :type fL: float
        :param fH: Upper cutoff frequency (Hz)
        :type fH: float
        :param order: Filter order
        :type order: int
        :param sm_wins: dictionary of smoothing windows in millisecond:
            - ppg: window for PPG signal
            - vpg: window for PPG' signal
            - apg: window for PPG" signal
            - jpg: window for PPG'" signal
        :type sm_wins: dict

        """

        self.fL = fL
        self.fH = fH
        self.order=order
        self.sm_wins=sm_wins

    def get_signals(self, s: DotMap):
        '''This function calculates the preprocessed PPG, PPG', PPG", and PPG'" signals.

        :param s: a struct of PPG signal:
            - s.v: a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.filtering: a bool for filtering
        :type s: DotMap

        :return: ppg, vpg, apg, jpg: preprocessed PPG, PPG', PPG", and PPG'"
        '''

        ## PPG filtering
        if s.filtering:
            fL = self.fL
            fH = self.fH
            order = self.order

            ppg_cb2 = _zero_phase_sos_filter(s.v, s.fs, fL, fH, order)

            if s.fs >= 75:
                win = _window_samples(s.fs, self.sm_wins['ppg'])
                ppg = _moving_average_zero_phase(ppg_cb2, win)
            else:
                ppg=ppg_cb2
        else:
            ppg=_as_finite_1d(s.v)

        if s.fs >= 150 and s.filtering:
            ## PPG' filtering
            win = _window_samples(s.fs, self.sm_wins['vpg'])
            dx = np.gradient(ppg)
            vpg = _moving_average_zero_phase(dx, win)

            ## PPG" filtering
            win = _window_samples(s.fs, self.sm_wins['apg'])
            ddx = np.gradient(vpg)
            apg = _moving_average_zero_phase(ddx, win)

            ## PPG'" filtering
            win = _window_samples(s.fs, self.sm_wins['jpg'])
            dddx = np.gradient(apg)
            jpg = _moving_average_zero_phase(dddx, win)
        else:
            vpg = np.gradient(ppg)
            apg = np.gradient(vpg)
            jpg = np.gradient(apg)

        return ppg, vpg, apg, jpg
