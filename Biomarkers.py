import pandas as pd
from scipy.signal import find_peaks
import scipy
import numpy as np

###########################################################################
############################ Get PPG Biomarkers ###########################
###########################################################################
def Biomarkers (sig, fs, fiducials):
    """
    This calc returns the main PPG biomarkers:
    CP, SUT, DT, SW50, DW50, DW50/SW50, Tpi, SA, SUT/CP, SOC, W50/Tpi, W50/SUT, SA/(Tpi-SUT), AUCPPG

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param fs: sampling frequency
    :type fs: int
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points
        PPG Fiducials Points.
        - Original signal: List of pulse onset, peak and dicrotic notch
        - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (a1,b1)
        - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a2, b2, c2, d2, e2)

    :return biomarkers: dictionary of biomarkers in different categories:
        - Original signal (Pulse onset, Pulse peak, Dicrotic notch)
        - Ratios of Systolic and Diastolic parts
        - 1st and 2nd derivative
        - Ratios of 1st and 2nd derivative’s points
    """

    BM_OSignal = get_BM_OSignal(sig, fs, fiducials)
    BM_ROSignal = get_BM_ROSignal(sig, fs, fiducials)
    BM_Derivatives = get_BM_Derivatives(sig, fs, fiducials)
    BM_RDerivatives = get_BM_RDerivatives(sig, fs, fiducials)

    biomarkers={'BM_OSignal': BM_OSignal , 'BM_ROSignal': BM_ROSignal, 'BM_Derivatives': BM_Derivatives, 'BM_RDerivatives': BM_RDerivatives}

    return biomarkers

###########################################################################
#################### Get Biomarkers of Original Signal ####################
###########################################################################
def get_BM_OSignal(sig, fs, fiducials):
    features_lst = ["CP",   # Cardiac Period, the time between two consecutive systolic peaks
                    "SUT",  # Systolic Upslope Time, the time between left onset and the systolic
                    "DT",   # Diastolic Time, the time between the systolic peak and right onset
                    "SW10", # Systolic Width, width at 10% of the pulse height from systolic part
                    "SW25", # Systolic Width, width at 25% of the pulse height from systolic part
                    "SW33", # Systolic Width, width at 33% of the pulse height from systolic part
                    "SW50", # Systolic Width, width at 50% of the pulse height from systolic part
                    "SW66", # Systolic Width, width at 66% of the pulse height from systolic part
                    "SW75", # Systolic Width, width at 75% of the pulse height from systolic part
                    "SW90", # Systolic Width, width at 90% of the pulse height from systolic part
                    "DW10", # Diastolic Width, width at 10% of the pulse height from diastolic part
                    "DW25", # Diastolic Width, width at 25% of the pulse height from diastolic part
                    "DW33", # Diastolic Width, width at 33% of the pulse height from diastolic part
                    "DW50", # Diastolic Width, width at 50% of the pulse height from diastolic part
                    "DW66", # Diastolic Width, width at 66% of the pulse height from diastolic part
                    "DW75", # Diastolic Width, width at 75% of the pulse height from diastolic part
                    "DW90", # Diastolic Width, width at 90% of the pulse height from diastolic part
                    "SW10+DW10", # Sum of Systolic and Diastolic Width at 10% width
                    "SW25+DW25", # Sum of Systolic and Diastolic Width at 25% width
                    "SW33+DW33", # Sum of Systolic and Diastolic Width at 33% width
                    "SW50+DW50", # Sum of Systolic and Diastolic Width at 50% width
                    "SW66+DW66", # Sum of Systolic and Diastolic Width at 66% width
                    "SW75+DW75", # Sum of Systolic and Diastolic Width at 75% width
                    "SW90+DW90", # Sum of Systolic and Diastolic Width at 90% width
                    "STT",       # Slope Transit Time, which based on geometrical considerations of  the PPG pulse wave to account for simultaneous
                    "AUCPPG",    # The area under the curve, a good indicator of change in vascular
                    "PIR",       # PPG Intensity Ratio, the ratio of Systolic Peak intensity and PPG valley intensity, reflects on the arterial diameter changes during one cardiac cycle from systole to diastole
                    "SA",        # Systolic Peak Amplitude
                    "SPT",       # Systolic Peak Time
                    "tpi",       # The time between the two onsets of the PPG systolic peak
                    "SOC"        # Systolic Peak Output Curve, Ratio between SUT and SA
                    #"DPT"       # Diastolic Peak Time
                    #"DNT"       # Dicrotic Notch Time
                    #"MS"        # Slope of the rising front is the amplitude of maximum upslope, normalized by the pulse amplitude
                    ]
    df, df_features = get_features(sig, fiducials['peaks'], fiducials['onsets'], fs, features_lst)

    return df_features

###########################################################################
################ Get Ratios of Systolic and Diastolic parts ###############
###########################################################################
def get_BM_ROSignal(sig, fs, fiducials):
    features_lst = ["DW10/SW10",    # Ratio of Systolic and Diastolic at 10% width
                    "DW25/SW25",    # Ratio of Systolic and Diastolic at 25% width
                    "DW33/SW33",    # Ratio of Systolic and Diastolic at 33% width
                    "DW50/SW50",    # Ratio of Systolic and Diastolic at 50% width
                    "DW66/SW66",    # Ratio of Systolic and Diastolic at 66% width
                    "DW75/SW75",    # Ratio of Systolic and Diastolic at 75% width
                    "DW90/SW90",    # Ratio of Systolic and Diastolic at 90% width
                    "SUT/CP",       # Ratio between SUT and CP
                    "SA/tpi-SUT",   # The ratio of SA and the difference between Tpi and SUT
                    "width_25_SUT", # The ratio of Width 25% and SUT
                    "width_25_tpi", # The ratio of Width 25% and Tpi
                    "width_50_SUT", # The ratio of Width 50% and SUT
                    "width_50_tpi", # The ratio of Width 50% and Tpi
                    "width_75_SUT", # The ratio of Width 75% and SUT
                    "width_75_tpi", # The ratio of Width 75% and Tpi
                    # "RI",           # Reflection Index is the ratio of DA and SA
                    # "SI",           # Stiffness Index is the ratio of SA and the difference between DPT and SUT
                    # "SC"            # Spring constant defined as x’’(sys) / (SA - MS) / SA, derived from a physical model of the elasticity of peripheral arteries
    ]

    df, df_features = get_features(sig, fiducials['peaks'], fiducials['onsets'], fs, features_lst)

    return df_features

###########################################################################
################# Get Biomarkers of 1st and 2nd Derivatives ###############
###########################################################################
def get_BM_Derivatives(sig, fs, fiducials):
    features_lst = ["a1",   # First maximum peak from 1st derivative of PPG waveform
                    "ta1",  # Time interval from the foot of PPG waveform to the time of with a1 occurs
                    "b1",   # First minimum peak from 1st derivative of PPG waveform
                    "tb1",  # Time interval from the foot of PPG waveform to the time of with b1 occurs
                    "a2",   # First maximum peak from 2nd derivative of PPG waveform
                    "ta2",  # Time interval from the foot of PPG waveform to the time of with a2 occurs
                    "b2",   # First minimum peak from 2nd derivative of PPG waveform
                    "tb2",  # Time interval from the foot of PPG waveform to the time of with b2 occurs
                    "c2",   # Second maximum peak from 2nd derivative of PPG waveform
                    "tc2",  # Time interval from the foot of PPG waveform to the time of with b2 occurs
                    "d2",   # Second minimum peak from 2nd derivative of PPG waveform
                    "td2",  # Time interval from the foot of PPG waveform to the time of with d1 occurs
                    "e2",   # Third maximum peak from 2nd derivative of PPG waveform
                    "te2",  # Time interval from the foot of PPG waveform to the time of with b1 occurs6
    ]

    df, df_features = get_features(sig, fiducials['peaks'], fiducials['onsets'], fs, features_lst)

    return df_features


###########################################################################
############### Get Ratios of 1st and 2nd derivative’s points #############
###########################################################################
def get_BM_RDerivatives(sig, fs, fiducials):
    features_lst = ["b1/a1",        # The ratio between minimum and maximum peaks of 1st PPG derivative
                    "ta1/cp",       # Ratio between ta1 and CP
                    "tb1/cp",       # Ratio between tb1 and CP
                    "ta2/cp",       # Ratio between ta2 and CP
                    "b2/a2",        # The ratio between minimum and maximum peaks of the second PPG derivative
                    "tb2/cp",       # Ratio between tb2 and CP
                    "c2/a2",        # The ratio between second maximum and maximum peaks of the second PPG derivative
                    "tc2/cp",       # Ratio between tc2 and CP
                    "d2/a2",        # The ratio between second minimum and maximum peaks of the second PPG derivative
                    "td2/cp",       # Ratio between td2 and CP
                    "e2/a2",        # The ratio between third maximum and maximum peaks of the second PPG derivative
                    "te2/cp",       # Ratio between te2 and CP
                    "(ta1-ta2)/cp", # The ratio between the interval maximum/minimum peaks of 1st derivative and CP
                    "(tb1-tb2)/cp", # The ratio between the interval
                    "(b2-c2-d2-e2)/a2", # Aging index of (b2-c2-d2-e2)/a2
                    "(b2-e2)=a2"        # Aging index of (b2-e2)/a2, instead of (b2-c2-d2-e2)/a2, when the c and d waves are missing
    ]

    df, df_features = get_features(sig, fiducials['peaks'], fiducials['onsets'], fs, features_lst)

    return df_features

###########################################################################
######################### PPG feature extraction ##########################
###########################################################################

class features_extract_PPG:

    def __init__(self, segment, peak_value, peak_time, next_peak_value, next_peak_time, onsets_values, onsets_times,
                 sample_rate,list_features):
        """
        :param segment: segment of PPG timeseries to analyse and extract features as a np array
        :param peak_value: PPG peak value
        :param peak_time: the time corresponding to the peak detected
        :param next_peak_value: PPG next peak value
        :param next_peak_time: the time corresponding to the peak detected
        :param onsets_values: array of PPG two onsets values surrounding the peak
        :param onsets_times: array of the two times corresponding to each onset detected
        :param sample_rate: segment data sample rate
        """

        self.list_features=list_features
        self.segment = segment
        self.peak_value = peak_value
        self.peak_time = peak_time
        self.next_peak_value = next_peak_value
        self.next_peak_time = next_peak_time
        self.onsets_values = onsets_values
        self.onsets_times = onsets_times
        self.sample_rate = sample_rate
        self.a1, self.b1, self.Ta1, self.Tb1 = self._getFirstDerivitivePoints()

        ## modified 11th of Nov. 2022
        self.a2, self.b2, self.c2, self.d2, self.e2, self.Ta2, self.Tb2, self.Tc2, self.Td2, self.Te2 = self._getSecondDerivitivePoints()
 #       self.a2, self.b2, self.Ta2, self.Tb2 = self._getSecondDerivitivePoints()

    def map_func(self):
        """ This function assign for each name of features a function that calculates it
            :returns my_funcs: a dictionary where the key is the name of the feature and
                               the value is the function to call"""
        my_funcs = {
                    "CP": self.getCP(),
                   "SUT": self.getSUT(),
                   "DT": self.getDT(),
                   "SW10": self.getSystolicWidth_d_percent(10),
                   "SW25": self.getSystolicWidth_d_percent(25),
                   "SW33": self.getSystolicWidth_d_percent(33),
                   "SW50": self.getSystolicWidth_d_percent(50),
                   "SW66": self.getSystolicWidth_d_percent(66),
                   "SW75": self.getSystolicWidth_d_percent(75),
                   "SW90": self.getSystolicWidth_d_percent(90),
                   "DW10": self.getDiastolicWidth_d_percent(10),
                   "DW25": self.getDiastolicWidth_d_percent(25),
                   "DW33": self.getDiastolicWidth_d_percent(33),
                   "DW50": self.getDiastolicWidth_d_percent(50),
                   "DW66": self.getDiastolicWidth_d_percent(66),
                   "DW75": self.getDiastolicWidth_d_percent(75),
                   "DW90": self.getDiastolicWidth_d_percent(90),
                    "DW10/SW10": self.getRatioSW_DW(10),
                    "DW25/SW25": self.getRatioSW_DW(25),
                    "DW33/SW33": self.getRatioSW_DW(33),
                    "DW50/SW50": self.getRatioSW_DW(50),
                    "DW66/SW66": self.getRatioSW_DW(66),
                    "DW75/SW75": self.getRatioSW_DW(75),
                    "DW90/SW90": self.getRatioSW_DW(90),
                    "SW10+DW10": self.getSumSW_DW(10),
                    "SW25+DW25": self.getSumSW_DW(25),
                    "SW33+DW33": self.getSumSW_DW(33),
                    "SW50+DW50": self.getSumSW_DW(50),
                    "SW66+DW66": self.getSumSW_DW(66),
                    "SW75+DW75": self.getSumSW_DW(75),
                    "SW90+DW90": self.getSumSW_DW(90),
                    "STT": self.getSTT(),
                    "AUCPPG": self.getAUCPPG(),
                    "PIR": self.getPIR(),
                    "SA": self.getSystolicPeak(),
                    "SPT": self.getSystolicPeakTime(),
                    "tpi": self.getTpi(),
                    "SOC": self.getSystolicPeakOutputCurve(),
                    "SUT/CP": self.getRatioSUTCP(),
                    "SA/tpi-SUT": self.getRatioSysPeakTpiSysTime(),
                    "width_25_SUT": self.getRatioWidth_SUT(25),
                    "width_25_tpi": self.getRatioWidth_tpi(25),
                    "width_50_SUT": self.getRatioWidth_SUT(50),
                    "width_50_tpi": self.getRatioWidth_tpi(50),
                    "width_75_SUT": self.getRatioWidth_SUT(75),
                    "width_75_tpi": self.getRatioWidth_tpi(75),
                    "a1": self.get_a1(),
                    "ta1": self.get_ta1(),
                    "b1": self.get_b1(),
                    "tb1": self.get_tb1(),
                    "a2": self.get_a2(),
                    "ta2": self.get_ta2(),
                    "b2": self.get_b2(),
                    "tb2": self.get_tb2(),
                    "c2": self.get_c2(),
                    "tc2": self.get_tc2(),
                    "d2": self.get_d2(),
                    "td2": self.get_td2(),
                    "e2": self.get_e2(),
                    "te2": self.get_te2(),
                    "ta1/cp": self.get_ratio_ta1_CP(),
                    "b1/a1": self.get_ratio_a1_b1(),
                    "tb1/cp": self.get_ratio_tb1_CP(),
                    "ta2/cp": self.get_ratio_ta2_CP(),
                    "b2/a2": self.get_ratio_a2_b2(),
                    "tb2/cp": self.get_ratio_tb2_CP(),
                    "c2/a2": self.get_ratio_a2_c2(),
                    "tc2/cp": self.get_ratio_tc2_CP(),
                    "d2/a2": self.get_ratio_a2_d2(),
                    "td2/cp": self.get_ratio_td2_CP(),
                    "e2/a2": self.get_ratio_a2_e2(),
                    "te2/cp": self.get_ratio_te2_CP(),
                    "(ta1-ta2)/cp": self.get_ratio_ta1_ta2_cp(),
                    "(tb1-tb2)/cp": self.get_ratio_tb1_tb2_cp(),
                    "(b2-c2-d2-e2)/a2": self.get_aging_index1(),
                    "(b2-e2)=a2": self.get_aging_index1()
        }
        return my_funcs

    def get_feature_extract_func(self):
        """ This function go through the list of features and call the function that is relevant to calculate it
            each feature takes two spots in the feature vector: one for the average value and one for its standard
            deviation.
            :returns features_vec: a vector of each feature avg value and std of the patient"""
        my_funcs = self.map_func()
        features_vec = []

        for feature in self.list_features:
            func_to_call = my_funcs[feature]
            features_vec.append(func_to_call)
        return features_vec

    def _getPeaksOnsets(self,x):
        """Find the peaks and onsets of a short FILTERED segment of PPG
        :return peaks
        :return onsets
        """
        peaks, _ = find_peaks(x)
        onsets, _ = find_peaks(-x)
        return peaks, onsets

    def _getFirstDerivitivePoints(self):
        """Calculate first derivitive points Ta1 and Tb1 from a SINGLE Onset-Onset segment of PPG
        :return Ta1: Time from PPG onset to peak of 1st derivitive
        :return Tb1: Time from PPG onset to onset of 1st derivitive
        """

        ## modified 11th of Nov. 2022
        sig = self.segment
        kernel_size = round(self.sample_rate / 20)
        kernel = np.ones(kernel_size) / kernel_size
        ma_sig = np.convolve(sig, kernel, mode='same')
        dx = np.gradient(ma_sig)

        #dx = np.gradient(self.segment)
        peaks, onsets = self._getPeaksOnsets(dx)
        a1 = peaks[0]
        b1 = onsets[0]
        pa1 = dx[a1]
        pb1 = dx[b1]
        Ta1 = peaks[0] / self.sample_rate
        Tb1 = onsets[0] / self.sample_rate

        return pa1, pb1, Ta1, Tb1

    def _getSecondDerivitivePoints(self):
        """Calculate second derivitive points Ta2 and Tb2 from a SINGLE Onset-Onset segment of PPG
        :return Ta1: Time from PPG onset to peak of 2st derivitive
        :return Tb1: Time from PPG onset to onset of 2st derivitive
        """

        ## modified 11th of Nov. 2022
        sig=self.segment
        kernel_size = round(self.sample_rate / 20)
        kernel = np.ones(kernel_size) / kernel_size
        ma_sig = np.convolve(sig, kernel, mode='same')
        dx = np.gradient(ma_sig)
        ma_dx = np.convolve(dx, kernel, mode='same')
        ddx = np.gradient(ma_dx)

        #dx = np.gradient(self.segment)
        #ddx = np.gradient(dx)
        peaks, onsets = self._getPeaksOnsets(ddx)
        a2 = peaks[0]
        b2 = onsets[0]
        pa2 = ddx[a2]
        pb2 = ddx[b2]
        Ta2 = peaks[0] / self.sample_rate
        Tb2 = onsets[0] / self.sample_rate

        ## modified 11th of Nov. 2022
        c2 = peaks[1]
        d2 = onsets[1]
        pc2 = ddx[c2]
        pd2 = ddx[d2]
        Tc2 = peaks[1] / self.sample_rate
        Td2 = onsets[1] / self.sample_rate

        if len(peaks)>2:
            e2 = peaks[2]
            pe2 = ddx[e2]
            Te2 = peaks[2] / self.sample_rate
        else:
            e2 = peaks[1]
            pe2 = ddx[e2]
            Te2 = peaks[1] / self.sample_rate

        return pa2, pb2, pc2, pd2, pe2,Ta2, Tb2, Tc2, Td2, Te2
        #return pa2, pb2, Ta2, Tb2
        ##

    def _find_nearest(self, arr, value):
        """ This function calculates the index in an array of the closest value to the arg value
            :param arr: the array where to find the index
            :param value: the value to be compared with
            :return idx: the index of the value closest to the arg value
            """
        idx = (np.abs(arr - value)).argmin()
        return idx

    def _getTime(self, vec, val):
        """ get the time of a value in the PPG waveform  data vector
            :param vec: the array where to find the index
            :param val: the value to be compared with
            :return index: the index of the value closest to the arg value
            """
        tmp_vec = np.array([vec[i]-val for i in range(0, len(vec))])
        index = self._find_nearest(tmp_vec, 0)
        return index

    def _getSysTime_from_val(self, val):
        """ get the time of a value in the PPG waveform  data vector
            :param val: the value from which we need the time in the timeserie
            :return t_data: the time corresponding to the value
            """
        idx_ons = int(self.onsets_times[0]*self.sample_rate)
        idx_peak = int(self.peak_time*self.sample_rate)
        idx = idx_peak - idx_ons
        vec_data = np.array(self.segment[0: idx])
        t_val = self._getTime(vec_data, val)
        t_data = t_val + idx_ons
        return t_data

    def _getDiaTime_from_val(self, val):
        """ get the time of a value in the PPG waveform  data vector
            :param val: the value from which we need the time in the timeserie
            :return t_data: the time corresponding to the value
            """
        idx_ons_right = int(self.onsets_times[1]*self.sample_rate)
        idx_ons_left = int(self.onsets_times[0] * self.sample_rate)
        idx_peak = int(self.peak_time*self.sample_rate)
        idx = idx_peak - idx_ons_left
        idx2 = idx_ons_right - idx_ons_left
        vec_data = np.array(self.segment[idx: idx2])
        t_val = self._getTime(vec_data, val)
        t_data = t_val + idx_peak
        return t_data

    def _getBaselineSlope(self):
        """ get the baseline slope in the PPG waveform  data vector
            :return baseline slope:
        """
        left_onset_time = self.onsets_times[0]
        right_onset_time = self.onsets_times[1]
        left_onset_value = self.onsets_values[0]
        right_onset_value = self.onsets_values[1]
        slope_numerator = right_onset_value - left_onset_value
        slope_denom = right_onset_time*self.sample_rate - left_onset_time*self.sample_rate
        return slope_numerator/slope_denom

    def _getBaselineCst(self):
        """ get the difference between the right and the left onset values
            :return cst:
        """
        left_onset_value = self.onsets_values[0]
        right_onset_value = self.onsets_values[1]
        cst = right_onset_value - left_onset_value
        return cst

    def getCP(self):
        """ CP means cardiac period and is the time difference between two peaks in a PPG waveform
            :return  CP feature:
        """
        cardiac_period = self.next_peak_time - self.peak_time
        return cardiac_period


    def getSUT(self):
        """ SUT means systolic upslope time and is the time difference between a peak and its left onset in a PPG waveform
        :return SUT feature:
        """
        left_onset = self.onsets_times[0]
        sut = self.peak_time - left_onset
        return sut

    def getDT(self):
        """ DT which means diastolic time is the time difference between a peak and its right onset in a PPG waveform
            :return DT feature:
        """
        right_onset = self.onsets_times[1]
        dt = right_onset - self.peak_time
        return dt


    def getSystolicWidth_d_percent(self, d):
        """ SW which means systolic width calculates the width of PPG waveform at d percent of the oulse height
            :return SW_d feature:
        """
        # value in segment corresponding to d percent of pulse height
        d_percent_val = (d/100)*(self.peak_value - self.onsets_values[0]) + self.onsets_values[0]
        time_of_d = self._getSysTime_from_val(d_percent_val)
        sw_d = self.peak_time - (time_of_d/self.sample_rate)
        return sw_d


    def getDiastolicWidth_d_percent(self, d):
        """ DW which means diastolic width calculates the width of PPG waveform at d percent of the pulse height
            :param d: the percentage chosen to calculate the width
            :return DW_d feature:
        """
        # value in segment corresponding to d percent of pulse height
        d_percent_val = (d/100)*(self.peak_value - self.onsets_values[1]) + self.onsets_values[1]
        time_of_d = self._getDiaTime_from_val(d_percent_val)
        dw_d = (time_of_d/self.sample_rate) - self.peak_time
        return dw_d

    def getSumSW_DW(self, d):
        """ The function calculates the sum of systolic and diastolic width at d percent of the pulse height
            :param d: the percentage chosen to calculate the width
            :return sum feature:
        """
        sw_d = self.getSystolicWidth_d_percent(d)
        dw_d = self.getDiastolicWidth_d_percent(d)
        return sw_d + dw_d


    def getRatioSW_DW(self, d):
        """ The function calculates the ratio of systolic and diastolic width at d percent of the pulse height
            :param d: the percentage chosen to calculate the width
            :return ratio feature:
        """
        sw_d = self.getSystolicWidth_d_percent(d)
        dw_d = self.getDiastolicWidth_d_percent(d)
        ratio = dw_d/sw_d
        return ratio

    def getPIR(self):
        """ The function calculates the ratio between the peak value and the right onset value
            :return pir feature:
        """
        pir = self.peak_value/self.onsets_values[1]
        return pir


    def getAUCPPG(self):
        """ The function calculates the area under the curve of a PPG waveform
            :return aucppg feature:
        """
        left_onset_time = self.onsets_times[0]*self.sample_rate
        right_onset_time = self.onsets_times[1]*self.sample_rate
        baseline_shift_slope = self._getBaselineSlope()
        baseline_cst = self._getBaselineCst()
        vec_value_between_ons = self.segment
        num_t = len(vec_value_between_ons)
        baseline = baseline_shift_slope*self.peak_time*self.sample_rate + baseline_cst
        sum = 0
        for t in range(0, num_t):
            sum += vec_value_between_ons[t] - baseline_shift_slope*((t+left_onset_time)) + baseline_cst
        AUCPPG_peak_sum_mod = 10*sum/((self.peak_value - baseline) * (right_onset_time - left_onset_time))
        # AUCPPG_peak_sum = sum
        return AUCPPG_peak_sum_mod

    def getUpslope(self):
        """ The function calculates Systolic Upslope between the left onset and the systolic peak.
            :return Systolic Upslope:
        """
        left_onset_time = self.onsets_times[0]*self.sample_rate
        left_onset_value = self.onsets_values[0]
        slope_numer = self.peak_value-left_onset_value
        slope_denom = self.peak_time*self.sample_rate - left_onset_time
        return slope_numer/slope_denom

    def getdiffVal(self):
        """ The function calculates the time between the left onset and the systolic peak.
            :return left onset and systolic peak time:
        """
        left_onset_value = self.onsets_values[0]
        diff = self.peak_value - left_onset_value
        return diff

    def getSTT(self):
        """ STT means slope transit time, which based on geometrical considerations of
            the PPG pulse wave to account for simultaneous.
            :return STT feature:
        """
        upslope = self.getUpslope()
        A = self.getdiffVal()
        return A/upslope

    def getSystolicPeak(self):
        """ The function calculates the Systolic Peak Amplitude.
            :return Systolic Peak Amplitude feature:
        """
        sys_peak = self.peak_value - self.onsets_values[0]
        return sys_peak

    def getSystolicPeakTime(self):
        """ Systolic Peak Time means the distance between the consecutive Systolic Peaks
             :return Systolic Peak Times:
         """
        return self.peak_time - self.onsets_times[0]

    def getSystolicPeakOutputCurve(self):
        """Peak time divided by systolic amplitude"""
        sys_peak_time = self.getSystolicPeakTime()
        sys_amplitude = self.getSystolicPeak()
        return sys_peak_time/sys_amplitude

    def getTpi(self):
        """ Tpi which means the time between the two onsets of the PPG systolic peak.
            :return Tpi feature:
        """
        return self.onsets_times[1] - self.onsets_times[0]


    def getRatioSUTCP(self):
        """ t1/cp is ratio between SUT and CP.
            :return t1/cp feature:
        """
        t1 = self.getSystolicPeakTime()
        cp = self.getCP()
        # print(cp)
        return t1/cp


    def getRatioSysPeakTpiSysTime(self):
        """ The function calculates the ratio of Systolic Peak Amplitude and the difference between Tpi and SUT.
            :return SA/(Tpi-SUT):
        """
        tpi = self.getTpi()
        t1 = self.getSystolicPeakTime()
        sys_peak = self.getSystolicPeak()
        return sys_peak/(tpi-t1)

    def getRatioWidth_tpi(self, d):
        """ The function calculates the ratio of Systolic+Diastolic width at d percent of the pulse height and Tpi.
            :param d: the percentage chosen to calculate the width
            :return Systolic+Diastolic width and Tpi ratio:
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        tpi = self.getTpi()
        return width/tpi

    def getRatioWidth_SUT(self, d):
        """ The function calculates the ratio of Systolic+Diastolic width at d percent of the pulse height and Systolic Peak Time.
            :param d: the percentage chosen to calculate the width
            :return Systolic+Diastolic width and Systolic Peak Time ratio:
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        t1 = self.getSystolicPeakTime()
        return width/t1

    def get_a1(self):
        """ The a1 means the first maximum peak from 1st derivative of PPG waveform.
            :return a1 feature:
        """
        return self.a1

    def get_ta1(self):
        """ The Ta1 means the time interval from the foot of PPG waveform to the time of with a1 occurs.
            :return Ta1 feature:
        """
        return self.Ta1

    def get_b1(self):
        """ The b1 means the first minimum peak from 1st derivative of PPG waveform.
            :return b1 feature:
        """
        return self.b1

    def get_tb1(self):
        """ The Tb1 means the time interval  from the foot of PPG waveform to the time of with b1 occurs.
            :return Tb1 feature:
        """
        return self.Tb1

    def get_a2(self):
        """ The a2 means the first maximum peak from 2nd derivative of PPG waveform.
            :return a2 feature:
        """
        return self.a2

    def get_ta2(self):
        """ Ta2 means the time interval from the foot of PPG waveform to the time of with a1 occurs.
            :return Ta2 feature:
        """
        return self.Ta2

    def get_b2(self):
        """ The b2 means the first minimum peak from 2nd derivative of PPG waveform.
            :return b2 feature:
        """
        return self.b2

    def get_tb2(self):
        """ Tb2 means the time interval from the foot of PPG waveform to the time of with b1 occurs.
            :return Tb2 feature:
        """
        return self.Tb2

    def get_c2(self):
        """ The c2 means the first minimum peak from 2nd derivative of PPG waveform.
            :return c2 feature:
        """
        return self.c2

    def get_tc2(self):
        """ Tc2 means the time interval from the foot of PPG waveform to the time of with b1 occurs.
            :return Tb2 feature:
        """
        return self.Tc2

    def get_d2(self):
        """ The d2 means the first minimum peak from 2nd derivative of PPG waveform.
            :return b2 feature:
        """
        return self.d2

    def get_td2(self):
        """ Td2 means the time interval from the foot of PPG waveform to the time of with b1 occurs.
            :return Td2 feature:
        """
        return self.Td2

    def get_e2(self):
        """ The e2 means the first minimum peak from 2nd derivative of PPG waveform.
            :return e2 feature:
        """
        return self.e2

    def get_te2(self):
        """ Te2 means the time interval from the foot of PPG waveform to the time of with b1 occurs.
            :return Te2 feature:
        """
        return self.Te2

    def get_ratio_ta1_CP(self):
        """ The function calculates the ratio of Ta1 and CP.
            :return Ta1 and CP ratio:
        """
        t1 = self.get_ta1()
        return t1 / self.getCP()

    def get_ratio_a1_b1(self):
        """ The a1 and b1 ratio is between minimum and maximum peaks of 1st PPG derivative.
            :return a1 and b1 ratio:
        """
        max_val = self.get_a1()
        min_val = self.get_b1()
        return min_val / max_val

    def get_ratio_tb1_CP(self):
        """ The function calculates the ratio of Tb1 and CP.
            :return Tb1 and CP ratio:
        """
        tb1 = self.get_tb1()
        return tb1 / self.getCP()

    def get_ratio_ta2_CP(self):
        """ The function calculates the ratio of ta2 and CP.
            :return ta2 and CP ratio feature:
        """
        ta2 = self.get_ta2()
        return ta2 / self.getCP()

    def get_ratio_a2_b2(self):
        """ The a2 and b2 ratio means ratio between minimum and maximum peaks of second PPG derivative.
            An important feature that reflects increased arterial stiffness. Related to distensibility of
            the peripheral artery and suggested that it is a useful index of atherosclerosis and altered
            arterial distensibility.
            :return a2 b2 ratio feature:
        """
        max_val = self.get_a2()
        min_val = self.get_b2()
        return min_val / max_val

    def get_ratio_tb2_CP(self):
        """ The function calculates the ratio of tb2 and CP.
            :return tb2 and CP ratio feature:
        """
        tb2 = self.get_tb2()
        return tb2 / self.getCP()

    def get_ratio_a2_c2(self):
        """ The ratio between second maximum and maximum peaks of the second PPG derivative.
            :return a2 c2 ratio feature:
        """
        max_val1 = self.get_a2()
        max_val2 = self.get_c2()
        return max_val1 / max_val2

    def get_ratio_tc2_CP(self):
        """ The function calculates the ratio of tc2 and CP.
            :return tc2 and CP ratio feature:
        """
        tc2 = self.get_tc2()
        return tc2 / self.getCP()

    def get_ratio_a2_d2(self):
        """ The ratio between second minimum and maximum peaks of the second PPG derivative.
            :return a2 d2 ratio feature:
        """
        max_val = self.get_a2()
        min_val = self.get_d2()
        return max_val / min_val

    def get_ratio_td2_CP(self):
        """ The function calculates the ratio of td2 and CP.
            :return td2 and CP ratio feature:
        """
        td2 = self.get_td2()
        return td2 / self.getCP()

    def get_ratio_a2_e2(self):
        """ The ratio between third maximum and maximum peaks of the second PPG derivative.
            :return a2 e2 ratio feature:
        """
        max_val1 = self.get_a2()
        max_val2 = self.get_e2()
        return max_val1 / max_val2

    def get_ratio_te2_CP(self):
        """ The function calculates the ratio of te2 and CP.
            :return te2 and CP ratio feature:
        """
        te2 = self.get_te2()
        return te2 / self.getCP()

    def get_ratio_ta1_ta2_cp(self):
        """ The function calculates the ratio between the interval maximum/minimum peaks of 1st derivative and CP
            :return (ta1 - ta2) / cp:
        """
        ta2 = self.get_ta2()
        ta1 = self.get_ta1()
        cp = self.getCP()
        return (ta1 - ta2) / cp

    def get_ratio_tb1_tb2_cp(self):
        """ The function calculates the ratio between the interval
            :return (tb1 - tb2) / cp:
        """
        tb2 = self.get_tb2()
        return (self.get_tb1() - tb2) / self.getCP()

    def get_aging_index1(self):
        """ The function calculates the Aging index of (b2-c2-d2-e2)/a2
            :return (b2-c2-d2-e2)/a2:
        """
        a2 = self.get_a2()
        b2 = self.get_b2()
        c2 = self.get_c2()
        d2 = self.get_d2()
        e2 = self.get_e2()
        return (b2-c2-d2-e2)/a2

    def get_aging_index2(self):
        """ The function calculates the Aging index of (b2-e2)/a2, instead of (b2-c2-d2-e2)/a2, when the c and d waves are missing
            :return (b2-e2)/a2:
        """
        a2 = self.get_a2()
        b2 = self.get_b2()
        e2 = self.get_e2()
        return (b2-e2)/a2

###########################################################################
############################# Get PPG features ############################
###########################################################################

def get_features(ppg, peaks, onsets, fs, features_lst):
    """
    The function calculates the biomedical features of PPG signal.

    :param ppg: 1-d array, of shape (N,) where N is the length of the signal
    :param peaks: 1-d array, peaks of the signal
    :param onsets: 1-d array, onsets of the signal
    :param fs: sampling frequency
    :type fs: int
    :param features_lst: list of features

    :return
        - df: data frame with onsets, offset and peaks
        - df_features: data frame with PPG signal features
    """

    df = pd.DataFrame()
    df_features = pd.DataFrame(columns=features_lst)
    # display(df_features)
    for i in range(len(onsets) - 1):
        #         #     print(f'i is {i}')
        onset = onsets[i]
        offset = onsets[i + 1]
        data = ppg[int(onset):int(offset)]
        peak = peaks[(peaks > onset) * (peaks < offset)]
        if len(peak) != 1:
            continue
        peak = peak[0]
        peak_value = ppg[peak]
        peak_time = peak / fs
        onset_value = ppg[onset]
        onset_time = onset / fs

        if (peak_value - onset_value) == 0:
            continue
        #     print(onset_time)
        offset_value = ppg[offset]
        offset_time = offset / fs
        #     print(offset_time)
        idx_array = np.where(peaks == peak)
        idx = idx_array[0]
        onsets_values = np.array([onset_value, offset_value])
        onsets_times = np.array([onset_time, offset_time])
        #     print(onset,peak, offset)
        if (idx + 1) < len(peaks):
            next_peak_value = ppg[peaks[idx + 1]][0]
            next_peak_time = peaks[idx + 1] / fs
            next_peak_time = next_peak_time[0]
            #         plt.plot(data)
            #         plt.show()
            #         print(peak_value,peak_time,next_peak_value,next_peak_time,onsets_values,onsets_times)
            try:
                features_extractor = features_extract_PPG(data, peak_value, peak_time, next_peak_value, next_peak_time,
                                                          onsets_values, onsets_times, fs, features_lst)
                features_vec = features_extractor.get_feature_extract_func()
                lst = list(features_vec)
                df_features.loc[len(df_features.index)] = lst
                #         display(df_features)
                df = df.append({'onset': onset, 'offset': offset, 'peak': peak}, ignore_index=True)
            except:
                pass
        else:
            print("no more peaks")
    return df, df_features