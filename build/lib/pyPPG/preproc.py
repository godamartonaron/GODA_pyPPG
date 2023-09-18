import numpy as np
from dotmap import DotMap
from scipy.signal import filtfilt
from scipy import signal

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
            - ppg: windows for PPG signal
            - vpg: windows for PPG' signal
            - apg: windows for PPG" signal
            - jpg: windows for PPG'" signal
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

            if fL==0:
                b,a = signal.cheby2(order, 20, [fH], 'low', fs=s.fs)
            else:
                b, a = signal.cheby2(order, 20, [fL,fH], 'bandpass', fs=s.fs)

            ppg_cb2 = filtfilt(b, a, s.v)

            if s.fs >= 75:
                win = round(s.fs * self.sm_wins['ppg']/1000)
                B = 1 / win * np.ones(win)
                ppg = filtfilt(B, 1, ppg_cb2)
            else:
                ppg=ppg_cb2
        else:
            ppg=s.v

        if s.fs >= 150 and s.filtering:
            ## PPG' filtering
            win = round(s.fs * self.sm_wins['vpg']/1000)
            B1 = 1 / win * np.ones(win)
            dx = np.gradient(ppg)
            vpg = filtfilt(B1, 1, dx)

            ## PPG" filtering
            win = round(s.fs * self.sm_wins['apg']/1000)
            B2 = 1 / win * np.ones(win)
            ddx = np.gradient(vpg)
            apg = filtfilt(B2, 1, ddx)

            ## PPG'" filtering
            win = round(s.fs * self.sm_wins['jpg']/1000)
            B3 = 1 / win * np.ones(win)
            dddx = np.gradient(apg)
            jpg = filtfilt(B3, 1, dddx)
        else:
            vpg = np.gradient(ppg)
            apg = np.gradient(vpg)
            jpg = np.gradient(apg)

        return ppg, vpg, apg, jpg