import copy
import pandas as pd

from pack_ppg._ErrorHandler import _check_shape_, WrongParameter
import numpy as np
from dotmap import DotMap
from scipy.signal import kaiserord, firwin, filtfilt, detrend, periodogram, lfilter, find_peaks, firls, resample
import matplotlib.pyplot as plt
import time
from scipy import signal

def Prefiltering(s):
    '''The function calculates the prefiltered PPG, PPg', PPG", and PPG'".

    :param s.v: a vector of PPG values
    :param s.fs: the sampling frequency of the PPG in Hz

    :return prefiltered PPG, PPg', PPG", and PPG'".
    '''

    ## PPG filtering
    fL = 0.5
    fH = 12
    order = 4
    b,a = signal.cheby2(order, 20, [fL, fH], 'bandpass', fs=s.fs)
    filt_sig_cb2 = filtfilt(b, a, s.v)

    if s.fs >= 75:
        win = round(s.fs * 0.02)
        B = 1 / win * np.ones(win)
        filt_sig = filtfilt(B, 1, filt_sig_cb2)
    else:
        filt_sig=filt_sig_cb2

    if s.fs >= 150:
        ## PPG' filtering
        win = round(s.fs * 0.01)
        B1 = 1 / win * np.ones(win)
        dx = np.gradient(s.v)
        filt_d1 = filtfilt(B1, 1, dx)

        ## PPG" filtering
        win = round(s.fs * 0.005)
        B2 = 1 / win * np.ones(win)
        dx = np.gradient(s.v)
        dx = filtfilt(B, 1, dx)
        ddx = np.gradient(dx)
        filt_d2 = filtfilt(B2, 1, ddx)

        ## PPG'" filtering
        filt_d3 = np.gradient(ddx)
    else:
        filt_d1 = np.gradient(filt_sig)
        filt_d2 = np.gradient(filt_d1)
        filt_d3 = np.gradient(filt_d2)

    return filt_sig, filt_d1, filt_d2, filt_d3