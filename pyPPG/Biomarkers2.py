import pandas as pd
from scipy.signal import find_peaks
import scipy
import numpy as np

###########################################################################
############################ Get PPG Biomarkers ###########################
###########################################################################
def Biomarkers2 (sig, fs, fiducials):
    """
    This calc returns the main PPG biomarkers:
    CP, SUT, DT, SW50, DW50, DW50/SW50, Tpi, SPA, SUT/CP, SOC, W50/Tpi, W50/SUT, SPA/(Tpi-SUT), AUCPPG

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param fs: sampling frequency
    :type fs: int
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points
        PPG Fiducials Points.
        - Original signal: List of pulse onset, peak and dicrotic notch
        - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)

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
                    "Tpi",  # The time between the two consecutive systolic onsets
                    "SUT",  # Systolic Upslope Time, the time between left systolic onset and the systolic peak
                    "Tsys", # The systolic time is between the systolic onsets and dicrotic notch
                    "Tdia", # The diastolic time is between the dicrotic notch and right onset
                    "DPT",  # Diastolic Peak Time, the time between the systolic onsets and diastolic peak
                    "dT",   # Time delay between systolic and diastolic peaks
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
                    "A1",        # Area bounded by times of pulse onset, between the systolic onsets and dicrotic notch
                    "A2",        # Area bounded by times of dicrotic notch, between the dicrotic notch and right onset
                    "PIR",       # PPG Intensity Ratio, the ratio of Systolic Peak intensity and PPG valley intensity, reflects on the arterial diameter changes during one cardiac cycle from systole to diastole
                    "SPA",       # Systolic Peak Amplitude
                    "DPA",       # Diastolic Peak Amplitude
                    "SOC"        # Systolic Peak Output Curve, Ratio between SUT and SPA
                    ]
    df, df_features = get_features(sig, fiducials, fs, features_lst)

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
                    "SPA/Tpi-SUT",   # The ratio of SPA and the difference between Tpi and SUT
                    "width_25_SUT", # The ratio of Width 25% and SUT
                    "width_25_Tpi", # The ratio of Width 25% and Tpi
                    "width_50_SUT", # The ratio of Width 50% and SUT
                    "width_50_Tpi", # The ratio of Width 50% and Tpi
                    "width_75_SUT", # The ratio of Width 75% and SUT
                    "width_75_Tpi" # The ratio of Width 75% and Tpi
                    # "RI",           # Reflection Index is the ratio of DPA and SPA
                    # "SI",           # Stiffness Index is the ratio of SPA and the difference between DPT and SUT
                    # "SC"            # Spring constant defined as x’’(sys) / (SPA - MS) / SPA, derived from a physical model of the elasticity of peripheral arteries
    ]

    df, df_features = get_features(sig, fiducials, fs, features_lst)

    return df_features

###########################################################################
################# Get Biomarkers of 1st and 2nd Derivatives ###############
###########################################################################
def get_BM_Derivatives(sig, fs, fiducials):
    features_lst = ["Tu",  # Time interval from the systolic onset to the time of with u occurs on PPG'.
                    "Tv",  # Time interval from the systolic onset to the time of with v occurs on PPG'.
                    "Tw",  # Time interval from the systolic onset to the time of with w occurs on PPG'.
                    "Ta",  # Time interval from the systolic onset to the time of with a occurs on PPG".
                    "Tb",  # Time interval from the systolic onset to the time of with b occurs on PPG".
                    "Tc",  # Time interval from the systolic onset to the time of with c occurs on PPG".
                    "Td",  # Time interval from the systolic onset to the time of with d occurs on PPG".
                    "Te",  # Time interval from the systolic onset to the time of with e occurs on PPG".
                    "Tf",  # Time interval from the systolic onset to the time of with f occurs on PPG".
    ]

    df, df_features = get_features(sig, fiducials, fs, features_lst)

    return df_features


###########################################################################
############### Get Ratios of 1st and 2nd derivative’s points #############
###########################################################################
def get_BM_RDerivatives(sig, fs, fiducials):
    features_lst = ["v/u",        # The ratio between minimum and maximum peaks of 1st PPG derivative
                    "Tu/CP",       # Ratio between Tu and CP
                    "Tv/CP",       # Ratio between Tv and CP
                    "Ta/CP",       # Ratio between Ta and CP
                    "b/a",        # The ratio between minimum and maximum peaks of the second PPG derivative
                    "Tb/CP",       # Ratio between Tb and CP
                    "c/a",        # The ratio between second maximum and maximum peaks of the second PPG derivative
                    "Tc/CP",       # Ratio between Tc and CP
                    "d/a",        # The ratio between second minimum and maximum peaks of the second PPG derivative
                    "Td/CP",       # Ratio between Td and CP
                    "e/a",        # The ratio between third maximum and maximum peaks of the second PPG derivative
                    "Te/CP",       # Ratio between Te and CP
                    "(Tu-Ta)/CP", # The ratio between the interval maximum/minimum peaks of 1st derivative and CP
                    "(Tv-Tb)/CP", # The ratio between the interval
                    "(b-c-d-e)/a", # Aging index of (b-c-d-e)/a
                    "(b-e)=a"        # Aging index of (b-e)/a, instead of (b-c-d-e)/a, when the c and d waves are missing
    ]

    df, df_features = get_features(sig, fiducials, fs, features_lst)

    return df_features

###########################################################################
######################### PPG feature extraction ##########################
###########################################################################

class features_extract_PPG:

    def __init__(self, segment, peak_value, peak_time, next_peak_value, next_peak_time, onsets_values, onsets_times,
                 sample_rate,list_features, fiducials):
        """
        :param segment: segment of PPG timeseries to analyse and extract features as a np array
        :param peak_value: PPG peak value
        :param peak_time: the time corresponding to the peak detected
        :param next_peak_value: PPG next peak value
        :param next_peak_time: the time corresponding to the peak detected
        :param onsets_values: array of PPG two onsets values surrounding the peak
        :param onsets_times: array of the two times corresponding to each onset detected
        :param sample_rate: segment data sample rate
        :param list_features: list of features
        :param fiducials: location of fiducial points of the given pulse wave
        """

        self.fiducials=fiducials
        self.list_features=list_features
        self.segment = segment
        self.peak_value = peak_value
        self.peak_time = peak_time
        self.next_peak_value = next_peak_value
        self.next_peak_time = next_peak_time
        self.onsets_values = onsets_values
        self.onsets_times = onsets_times
        self.sample_rate = sample_rate
        self.dn, self.dp, self.Tdn, self.Tdp = self._getDicroticNotchDiastolicPeak()
        self.u, self.v, self.w, self.Tu, self.Tv, self.Tw = self._getFirstDerivitivePoints()
        self.a, self.b, self.c, self.d, self.e, self.f, self.Ta, self.Tb, self.Tc, self.Td, self.Te, self.Tf = self._getSecondDerivitivePoints()

    def map_func(self):
        """ This function assign for each name of features a function that calculates it
            :returns my_funcs: a dictionary where the key is the name of the feature and
                               the value is the function to call"""
        my_funcs = {"CP": self.getCP(),
                   "Tpi": self.getTpi(),
                   "SUT": self.getSUT(),
                   "Tsys": self.getTsys(),
                   "Tdia": self.getTdia(),
                   "DPT": self.getDPT(),
                   "dT": self.getdT(),
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
                    "A1": self.getA1(),
                    "A2": self.getA2(),
                    "PIR": self.getPIR(),
                    "SPA": self.getSystolicPeak(),
                    "DPA": self.getDiastolicPeak(),
                    "SPT": self.getSystolicPeakTime(),
                    "SOC": self.getSystolicPeakOutputCurve(),
                    "SUT/CP": self.getRatioSUTCP(),
                    "SPA/Tpi-SUT": self.getRatioSysPeakTpiSysTime(),
                    "width_25_SUT": self.getRatioWidth_SUT(25),
                    "width_25_Tpi": self.getRatioWidth_Tpi(25),
                    "width_50_SUT": self.getRatioWidth_SUT(50),
                    "width_50_Tpi": self.getRatioWidth_Tpi(50),
                    "width_75_SUT": self.getRatioWidth_SUT(75),
                    "width_75_Tpi": self.getRatioWidth_Tpi(75),
                    "u": self.get_u(),
                    "Tu": self.get_Tu(),
                    "v": self.get_v(),
                    "Tv": self.get_Tv(),
                    "w": self.get_v(),
                    "Tw": self.get_Tv(),
                    "a": self.get_a(),
                    "Ta": self.get_Ta(),
                    "b": self.get_b(),
                    "Tb": self.get_Tb(),
                    "c": self.get_c(),
                    "Tc": self.get_Tc(),
                    "d": self.get_d(),
                    "Td": self.get_Td(),
                    "e": self.get_e(),
                    "Te": self.get_Te(),
                    "f": self.get_e(),
                    "Tf": self.get_Te(),
                    "Tu/CP": self.get_ratio_Tu_CP(),
                    "v/u": self.get_ratio_u_v(),
                    "Tv/CP": self.get_ratio_Tv_CP(),
                    "Ta/CP": self.get_ratio_Ta_CP(),
                    "b/a": self.get_ratio_a_b(),
                    "Tb/CP": self.get_ratio_Tb_CP(),
                    "c/a": self.get_ratio_a_c(),
                    "Tc/CP": self.get_ratio_Tc_CP(),
                    "d/a": self.get_ratio_a_d(),
                    "Td/CP": self.get_ratio_Td_CP(),
                    "e/a": self.get_ratio_a_e(),
                    "Te/CP": self.get_ratio_Te_CP(),
                    "(Tu-Ta)/CP": self.get_ratio_Tu_Ta_CP(),
                    "(Tv-Tb)/CP": self.get_ratio_Tv_Tb_CP(),
                    "(b-c-d-e)/a": self.get_aging_index1(),
                    "(b-e)=a": self.get_aging_index1()
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

    def _getDicroticNotchDiastolicPeak(self):
        """Calculate Dicrotic Notch and Diastolic Peak of PPG
        :return dn: Sample distance from PPG onset to the Dicrotic Notch on PPG
        :return dp: Sample distance from PPG onset to the Diastolic Peak on PPG
        :return Tdn: Time from PPG onset to the Dicrotic Notch on PPG
        :return Tdp: Time from PPG onset to the Diastolic Peak on PPG
        """

        dn = (self.fiducials.dn-self.fiducials.os).values[0]
        dp = (self.fiducials.dp-self.fiducials.os).values[0]
        Tdn = dn / self.sample_rate
        Tdp = dp / self.sample_rate

        return dn, dp, Tdn, Tdp
    def _getFirstDerivitivePoints(self):
        """Calculate first derivitive points from a SINGLE Onset-Onset segment of PPG'
        :return u: Sample distance from PPG onset to the greatest maximum peak between the left systolic onset and the systolic peak on PPG'
        :return v: Sample distance from PPG onset to the lowest minimum pits between the systolic peak and the right systolic onset on PPG'
        :return w: Sample distance from PPG onset to the first maximum peak after v on PPG'
        :return Tu: Time from PPG onset to the greatest maximum peak between the left systolic onset and the systolic peak on PPG'
        :return Tv: Time from PPG onset to the lowest minimum pits between the systolic peak and the right systolic onset on PPG'
        :return Tw: Time from PPG onset to the first maximum peak after v on PPG'
        """

        u = (self.fiducials.u-self.fiducials.os).values[0]
        v = (self.fiducials.v-self.fiducials.os).values[0]
        w = (self.fiducials.w - self.fiducials.os).values[0]
        Tu = u / self.sample_rate
        Tv = v / self.sample_rate
        Tw = w / self.sample_rate

        return u, v, w, Tu, Tv, Tw

    def _getSecondDerivitivePoints(self):
        """Calculate second derivitive points from a SINGLE Onset-Onset segment of PPG"
        :return a: Sample distance from PPG onset to the first maximum peak between left systolic onset and systolic peak on PPG"
        :return b: Sample distance from PPG onset to the first minimum pits after a on PPG"
        :return c: Sample distance from PPG onset to the greatest maximum peak between b and e, or if no maximum peak is present then the inflection point on PPG"
        :return d: Sample distance from PPG onset to the lowest minimum pits between c and e, or if no minimum pits is present then the inflection point on PPG"
        :return e: Sample distance from PPG onset to the greatest maximum peak between the systolic peak and  the right systolic onset on PPG"
        :return f: Sample distance from PPG onset to the first minimum pits after e on PPG"
        :return Ta: Time from PPG onset to the first maximum peak between left systolic onset and systolic peak on PPG"
        :return Tb: Time from PPG onset to the first minimum pits after a on PPG"
        :return Tc: Time from PPG onset to the greatest maximum peak between b and e, or if no maximum peak is present then the inflection point on PPG"
        :return Td: Time from PPG onset to the lowest minimum pits between c and e, or if no minimum pits is present then the inflection point on PPG"
        :return Te: Time from PPG onset to the greatest maximum peak between the systolic peak and  the right systolic onset on PPG"
        :return Tf: Time from PPG onset to the first minimum pits after e on PPG"
        """

        a = (self.fiducials.a-self.fiducials.os).values[0]
        b = (self.fiducials.b-self.fiducials.os).values[0]
        c = (self.fiducials.c - self.fiducials.os).values[0]
        d = (self.fiducials.d-self.fiducials.os).values[0]
        e = (self.fiducials.e-self.fiducials.os).values[0]
        f = (self.fiducials.f - self.fiducials.os).values[0]
        Ta = a / self.sample_rate
        Tb = b / self.sample_rate
        Tc = c / self.sample_rate
        Td = d / self.sample_rate
        Te = e / self.sample_rate
        Tf = f / self.sample_rate

        return a, b, c, d, e, f, Ta, Tb, Tc, Td, Te, Tf

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

    def getTsys(self):
        """ Tsys means systolic the systolic time between the left onsets and dicrotic notch
        :return Tsys feature:
        """

        Tsys = self.Tdn
        return Tsys

    def getTdia(self):
        """ Tdia means systolic the diastolic time between the dicrotic notch and right onset
        :return Tdia feature:
        """

        Tdia = self.getTpi()-self.Tdp
        return Tdia

    def getDPT(self):
        """ DPT means the time between the left onsets and diastolic peak
        :return DPT feature:
        """

        DPT = self.Tdp
        return DPT

    def getdT(self):
        """ dT means time delay between systolic and diastolic peaks
        :return dT feature:
        """

        dT = self.Tdp-self.getSUT()
        return dT

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

    def getA1(self):
        """ The function calculates the area under the curve between the systolic onsets and dicrotic notch
            :return A1 feature:
        """
        left_onset_time = self.onsets_times[0]*self.sample_rate
        right_onset_time = self.onsets_times[1]*self.sample_rate
        baseline_shift_slope = self._getBaselineSlope()
        baseline_cst = self._getBaselineCst()
        vec_value_between_ons = self.segment
        num_t = self.dn
        baseline = baseline_shift_slope*self.peak_time*self.sample_rate + baseline_cst
        sum = 0
        for t in range(0, num_t):
            sum += vec_value_between_ons[t] - baseline_shift_slope*((t+left_onset_time)) + baseline_cst
        A1= 10*sum/((self.peak_value - baseline) * (right_onset_time - left_onset_time))

        return A1

    def getA2(self):
        """ The function calculates the area under the curve between the dicrotic notch and right onset
            :return A2 feature:
        """
        A2=self.getAUCPPG()-self.getA1()
        return A2

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

    def getDiastolicPeak(self):
        """ The function calculates the Diastolic Peak Amplitude.
            :return Diastolic Peak Amplitude feature:
        """
        ## temp solution 04/04/2023 --> if DN==DP
        dp_value = self.segment[self.dp+20]
        dia_peak = dp_value - self.onsets_values[0]
        return dia_peak

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
        """ T1/CP is ratio between SUT and CP.
            :return T1/CP feature:
        """
        T1 = self.getSystolicPeakTime()
        CP = self.getCP()
        # print(CP)
        return T1/CP


    def getRatioSysPeakTpiSysTime(self):
        """ The function calculates the ratio of Systolic Peak Amplitude and the difference between Tpi and SUT.
            :return SPA/(Tpi-SUT):
        """
        Tpi = self.getTpi()
        T1 = self.getSystolicPeakTime()
        sys_peak = self.getSystolicPeak()
        return sys_peak/(Tpi-T1)

    def getRatioWidth_Tpi(self, d):
        """ The function calculates the ratio of Systolic+Diastolic width at d percent of the pulse height and Tpi.
            :param d: the percentage chosen to calculate the width
            :return Systolic+Diastolic width and Tpi ratio:
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        Tpi = self.getTpi()
        return width/Tpi

    def getRatioWidth_SUT(self, d):
        """ The function calculates the ratio of Systolic+Diastolic width at d percent of the pulse height and Systolic Peak Time.
            :param d: the percentage chosen to calculate the width
            :return Systolic+Diastolic width and Systolic Peak Time ratio:
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        T1 = self.getSystolicPeakTime()
        return width/T1

    def get_u(self):
        """ The u means the sample interval from the systolic onset to the time of with u occurs on PPG'.
            :return u feature:
        """
        return self.u

    def get_Tu(self):
        """ The Tu means the time interval from the systolic onset to the time of with u occurs on PPG'.
            :return Tu feature
        """
        return self.Tu

    def get_v(self):
        """ The v means the sample interval from the systolic onset to the time of with v occurs on PPG'.
            :return v feature:
        """
        return self.v

    def get_Tv(self):
        """ The Tv means the time interval from the systolic onset to the time of with v occurs on PPG'.
            :return Tv feature:
        """
        return self.Tv

    def get_w(self):
        """ The w means the sample interval from the systolic onset to the time of with w occurs on PPG'.
            :return v feature:
        """
        return self.w

    def get_Tw(self):
        """ The Tv means the time interval from the systolic onset to the time of with w occurs on PPG'.
            :return Tw feature:
        """
        return self.Tw

    def get_a(self):
        """ The a means the sample interval from the systolic onset to the time of with a occurs on PPG".
            :return a feature:
        """
        return self.a

    def get_Ta(self):
        """ Ta means the time interval from the systolic onset to the time of with a occurs on PPG".
            :return Ta feature:
        """
        return self.Ta

    def get_b(self):
        """ The b means the sample interval from the systolic onset to the time of with b occurs on PPG".
            :return b feature:
        """
        return self.b

    def get_Tb(self):
        """ Tb means the time interval from the systolic onset to the time of with b occurs on PPG".
            :return Tb feature:
        """
        return self.Tb

    def get_c(self):
        """ The c means the sample interval from the systolic onset to the time of with c occurs on PPG".
            :return c feature:
        """
        return self.c

    def get_Tc(self):
        """ Tc means the time interval from the systolic onset to the time of with c occurs on PPG".
            :return Tb feature:
        """
        return self.Tc

    def get_d(self):
        """ The d means the sample interval from the systolic onset to the time of with d occurs on PPG".
            :return b feature:
        """
        return self.d

    def get_Td(self):
        """ Td means the time interval from the systolic onset to the time of with d occurs on PPG".
            :return Td feature:
        """
        return self.Td

    def get_e(self):
        """ The e means the sample interval from the systolic onset to the time of with e occurs on PPG".
            :return e feature:
        """
        return self.e

    def get_Te(self):
        """ Te means the time interval from the systolic onset to the time of with e occurs on PPG".
            :return Te feature:
        """
        return self.Te

    def get_f(self):
        """ The f means the sample interval from the systolic onset to the time of with f occurs on PPG".
            :return f feature:
        """
        return self.f

    def get_Tf(self):
        """ Tf means the time interval from the systolic onset to the time of with f occurs on PPG".
            :return Tf feature:
        """
        return self.Tf

    def get_ratio_Tu_CP(self):
        """ The function calculates the ratio of Tu and CP.
            :return Tu and CP ratio:
        """
        T1 = self.get_Tu()
        return T1 / self.getCP()

    def get_ratio_u_v(self):
        """ The u and v ratio is between minimum and maximum peaks of 1st PPG derivative.
            :return u and v ratio:
        """
        max_val = self.get_u()
        min_val = self.get_v()
        return min_val / max_val

    def get_ratio_Tv_CP(self):
        """ The function calculates the ratio of Tv and CP.
            :return Tv and CP ratio:
        """
        Tv = self.get_Tv()
        return Tv / self.getCP()

    def get_ratio_Ta_CP(self):
        """ The function calculates the ratio of Ta and CP.
            :return Ta and CP ratio feature:
        """
        Ta = self.get_Ta()
        return Ta / self.getCP()

    def get_ratio_a_b(self):
        """ The a and b ratio means ratio between minimum and maximum peaks of second PPG derivative.
            An important feature that reflects increased arterial stiffness. Related to distensibility of
            the peripheral artery and suggested that it is a useful index of atherosclerosis and altered
            arterial distensibility.
            :return a b ratio feature:
        """
        max_val = self.get_a()
        min_val = self.get_b()
        return min_val / max_val

    def get_ratio_Tb_CP(self):
        """ The function calculates the ratio of Tb and CP.
            :return Tb and CP ratio feature:
        """
        Tb = self.get_Tb()
        return Tb / self.getCP()

    def get_ratio_a_c(self):
        """ The ratio between second maximum and maximum peaks of the second PPG derivative.
            :return a c ratio feature:
        """
        max_val1 = self.get_a()
        max_val2 = self.get_c()
        return max_val1 / max_val2

    def get_ratio_Tc_CP(self):
        """ The function calculates the ratio of Tc and CP.
            :return Tc and CP ratio feature:
        """
        Tc = self.get_Tc()
        return Tc / self.getCP()

    def get_ratio_a_d(self):
        """ The ratio between second minimum and maximum peaks of the second PPG derivative.
            :return a d ratio feature:
        """
        max_val = self.get_a()
        min_val = self.get_d()
        return max_val / min_val

    def get_ratio_Td_CP(self):
        """ The function calculates the ratio of Td and CP.
            :return Td and CP ratio feature:
        """
        Td = self.get_Td()
        return Td / self.getCP()

    def get_ratio_a_e(self):
        """ The ratio between third maximum and maximum peaks of the second PPG derivative.
            :return a e ratio feature:
        """
        max_val1 = self.get_a()
        max_val2 = self.get_e()
        return max_val1 / max_val2

    def get_ratio_Te_CP(self):
        """ The function calculates the ratio of Te and CP.
            :return Te and CP ratio feature:
        """
        Te = self.get_Te()
        return Te / self.getCP()

    def get_ratio_Tu_Ta_CP(self):
        """ The function calculates the ratio between the interval maximum/minimum peaks of 1st derivative and CP
            :return (Tu - Ta) / CP:
        """
        Ta = self.get_Ta()
        Tu = self.get_Tu()
        CP = self.getCP()
        return (Tu - Ta) / CP

    def get_ratio_Tv_Tb_CP(self):
        """ The function calculates the ratio between the interval
            :return (Tv - Tb) / CP:
        """
        Tb = self.get_Tb()
        return (self.get_Tv() - Tb) / self.getCP()

    def get_aging_index1(self):
        """ The function calculates the Aging index of (b-c-d-e)/a
            :return (b-c-d-e)/a:
        """
        a = self.get_a()
        b = self.get_b()
        c = self.get_c()
        d = self.get_d()
        e = self.get_e()
        return (b-c-d-e)/a

    def get_aging_index2(self):
        """ The function calculates the Aging index of (b-e)/a, instead of (b-c-d-e)/a, when the c and d waves are missing
            :return (b-e)/a:
        """
        a = self.get_a()
        b = self.get_b()
        e = self.get_e()
        return (b-e)/a

###########################################################################
############################# Get PPG features ############################
###########################################################################

def get_features(ppg, fiducials, fs, features_lst):
    """
    The function calculates the biomedical features of PPG signal.

    :param ppg: 1-d array, of shape (N,) where N is the length of the signal
    :param fiducials: M-d Dateframe, where M is the number of fiducial points
    :param fs: sampling frequency
    :type fs: int
    :param features_lst: list of features

    :return
        - df: data frame with onsets, offset and peaks
        - df_features: data frame with PPG signal features
    """

    df = pd.DataFrame()
    df_features = pd.DataFrame(columns=features_lst)
    peaks = fiducials.pk.values
    onsets = fiducials.os.values

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

        temp_fiducials = fiducials.iloc[[i]]

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
            next_peak_value = ppg[peaks[idx + 1].astype('int64')][0]
            next_peak_time = peaks[idx + 1] / fs
            next_peak_time = next_peak_time[0]
            #         plt.plot(data)
            #         plt.show()
            #         print(peak_value,peak_time,next_peak_value,next_peak_time,onsets_values,onsets_times)
            try:
                features_extractor = features_extract_PPG(data, peak_value, peak_time, next_peak_value, next_peak_time,
                                                          onsets_values, onsets_times, fs, features_lst,temp_fiducials)
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