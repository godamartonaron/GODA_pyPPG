import pandas as pd
from scipy.signal import find_peaks
import scipy
import numpy as np
from dotmap import DotMap

###########################################################################
############################ Get PPG Biomarkers ###########################
###########################################################################
def Biomarkers2 (s, fiducials):
    """
    This calc returns the main PPG biomarkers:
    CP, SUT, DT, SW50, DW50, DW50/SW50, Tpi, SPA, SUT/CP, SOC, W50/Tpi, W50/SUT, SPA/(Tpi-SUT), AUCPPG

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points
        PPG Fiducials Points.
        - Original signal: List of pulse onset, peak and dicrotic notch
        - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
        - 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)

    :return biomarkers: dictionary of biomarkers in different categories:
        - Original signal (Pulse onset, Pulse peak, Dicrotic notch)
        - Ratios of Systolic and Diastolic parts
        - 1st and 2nd derivative
        - Ratios of 1st and 2nd derivative’s points
    """

    BM_OSignal = get_BM_OSignal(s, fiducials)
    BM_ROSignal = get_BM_ROSignal(s, fiducials)
    BM_Derivatives = get_BM_Derivatives(s, fiducials)
    BM_RDerivatives = get_BM_RDerivatives(s, fiducials)

    biomarkers={'BM_OSignal': BM_OSignal , 'BM_ROSignal': BM_ROSignal, 'BM_Derivatives': BM_Derivatives, 'BM_RDerivatives': BM_RDerivatives}

    return biomarkers

###########################################################################
#################### Get Biomarkers of Original Signal ####################
###########################################################################
def get_BM_OSignal(s, fiducials):
    features_lst = ["Tpi",   # Pulse Interval, the time between the pulse onset and pulse offset
                    "Tpp",   # Peak-to-Peak Interval, the time between two consecutive systolic peaks
                    "Tsys",	 # Systolic Time, the time between the pulse onset and dicrotic notch
                    "Tdia",  # Diastolic Time, the time is between the dicrotic notch and pulse offset
                    "Tsp",   # Systolic Peak Time, the time between the pulse onset and systolic peak
                    "Tdp",	 # Diastolic Peak Time, the time between the pulse onset and diastolic peak
                    "deltaT",# Time Delay, the time between the systolic peak and diastolic peak
                    "Tsw10", # Systolic Width, the width at 10% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw25", # Systolic Width, the width at 25% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw33", # Systolic Width, the width at 33% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw50", # Systolic Width, the width at 50% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw66", # Systolic Width, the width at 66% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw75", # Systolic Width, the width at 75% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw90", # Systolic Width, the width at 90% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tdw10", # Diastolic Width, the width at 10% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw25", # Diastolic Width, the width at 25% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw33", # Diastolic Width, the width at 33% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw50", # Diastolic Width, the width at 50% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw66", # Diastolic Width, the width at 66% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw75", # Diastolic Width, the width at 75% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw90", # Diastolic Width, the width at 90% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tpw10", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 10%
                    "Tpw25", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 25%
                    "Tpw33", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 33%
                    "Tpw50", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 50%
                    "Tpw66", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 66%
                    "Tpw75", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 75%
                    "Tpw90", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 90%
                    "Asp",   # Systolic Peak Amplitude, the difference in amplitude between the pulse onset and systolic peak
                    "Adn",   # Dicrotic Notch Amplitude, the difference in amplitude between the pulse onset and dicrotic notch
                    "Adp",   # Diastolic Peak Amplitude, the difference in amplitude between the pulse onset and diastolic peak
                    "Aoff",  # Pulse Onset Amplitude, the difference in amplitude between the pulse onset and pulse offset
                    "AUCpi", # Area Under Pulse Interval Curve, the area under the pulse wave between pulse onset and pulse offset
                    "AUCsys",# Area Under Systolic Curve, the area under the pulse wave between the pulse onset and the dicrotic notch
                    "AUCdia",# Area Under Diastolic Curve, the area under the pulse wave between the dicrotic notch and pulse offset
                    ]
    df, df_features = get_features(s, fiducials, features_lst)

    return df_features

###########################################################################
################ Get Ratios of Systolic and Diastolic parts ###############
###########################################################################
def get_BM_ROSignal(s, fiducials):
    features_lst = ["IPR",          # Instantaneous Pulse Rate, 60 / Tpi
                    "Tsys/Tdia",    # The ratio of the Systolic Time to the Diastolic Time
                    "Tpw25/Tpi",    # The ratio of the Pulse Width at 25% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw50/Tpi",    # The ratio of the Pulse Width at 50% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw75/Tpi",    # The ratio of the Pulse Width at 75% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw25/Tsp",    # The ratio of the Pulse Width at 25% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tpw50/Tsp",    # The ratio of the Pulse Width at 50% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tpw75/Tsp",    # The ratio of the Pulse Width at 75% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tdw10/Tsw10",  # The ratio of the Diastolic Width to the Systolic Width at 10% width
                    "Tdw25/Tsw25",  # The ratio of the Diastolic Width to the Systolic Width at 25% width
                    "Tdw33/Tsw33",  # The ratio of the Diastolic Width to the Systolic Width at 33% width
                    "Tdw50/Tsw50",  # The ratio of the Diastolic Width to the Systolic Width at 50% width
                    "Tdw66/Tsw66",  # The ratio of the Diastolic Width to the Systolic Width at 66% width
                    "Tdw75/Tsw75",  # The ratio of the Diastolic Width to the Systolic Width at 75% width
                    "Tdw90/Tsw90",  # The ratio of the Diastolic Width to the Systolic Width at 90% width
                    "Tsp/Tpi",      # The ratio of the Systolic Peak Time to the Pulse Interval
                    "Asp/Aoff",     # The ratio of the Systolic Peak Amplitude to the Pulse Offset Amplitude
                    "Adp/Asp",      # Reflection Index, the ratio of the Diastolic Peak Amplitude to the Systolic Peak Amplitude
                    "IPA",          # Inflection Point Area, the ratio of the Area Under Diastolic Curve to the Area Under Systolic Curve
                    "Tsp/Asp",      # The ratio of the Systolic Peak Time to the Systolic Peak Amplitude
                    "Asp/deltaT",   # Stiffness Index, the ratio of the Systolic Peak Amplitude to the Time Delay
                    "Asp/(Tpi-Tsp)",# The ratio of the Systolic Peak Amplitude to the difference between the Pulse Interval and Systolic Peak Time
    ]

    df, df_features = get_features(s, fiducials, features_lst)

    return df_features

###########################################################################
################# Get Biomarkers of 1st and 2nd Derivatives ###############
###########################################################################
def get_BM_Derivatives(s, fiducials):
    features_lst = ["Tu",       # u-point time, the time between the pulse onset and u-point
                    "Tv",       # v-point time, the time between the pulse onset and v-point
                    "Tw",       # w-point time, the time between the pulse onset and w-point
                    "Ta",       # a-point time, the time between the pulse onset and a-point
                    "Tb",       # b-point time, the time between the pulse onset and b-point
                    "Tc",       # c-point time, the time between the pulse onset and c-point
                    "Td",       # d-point time, the time between the pulse onset and d-point
                    "Te",       # e-point time, the time between the pulse onset and e-point
                    "Tf",       # f-point time, the time between the pulse onset and f-point
                    "Tb–c",	    # b–c interval time, the time between the b-point and c-point
                    "Tb–d",	    # b–d interval time, the time between the b-point and d-point
                    "Tp1",	    # p1-point time, the time between the pulse onset and p1-point
                    "Tp2",      # p2-point time, the time between the pulse onset and p2-point
                    "Tp1–dp",   # p1–dia interval time, the time between the p1-point and diastolic peak
                    "Tp2–dp",   # p2–dia interval time, the time between the p2-point and diastolic peak
    ]

    df, df_features = get_features(s, fiducials, features_lst)

    return df_features


###########################################################################
############### Get Ratios of 1st and 2nd derivative’s points #############
###########################################################################
def get_BM_RDerivatives(s, fiducials):
    features_lst = ["Tu/CP",      # Ratio between Tu and CP
                    "Tv/CP",      # Ratio between Tv and CP
                    "Tw/CP",      # Ratio between Tw and CP
                    "Ta/CP",      # Ratio between Ta and CP
                    "Tb/CP",      # Ratio between Tb and CP
                    "Tc/CP",      # Ratio between Tc and CP
                    "Td/CP",      # Ratio between Td and CP
                    "Te/CP",      # Ratio between Te and CP
                    "Tf/CP",      # Ratio between Te and CP
                    "(Tu-Ta)/CP", # The ratio between the interval maximum/minimum peaks of 1st derivative and CP
                    "(Tv-Tb)/CP", # The ratio between the interval
                    "AGI_bcdef/a",# Aging index of (b-c-d-e-f)/a
                    "AGI_bcde/a", # Aging index of (b-c-d-e)/a
                    "AGI_bcd/a",  # Aging index of (b-c-d)/a
                    "AGI_be/a",   # Aging index of (b-e)/a, instead of (b-c-d-e)/a, when the c and d waves are missing
                    "v/u",  # The ratio v and u of the PPG'
                    "w/u",  # The ratio w and u of the PPG'
                    "b/a",  # The ratio between b and a of the PPG"
                    "c/a",  # The ratio between c and a of the PPG"
                    "d/a",  # The ratio between d and a of the PPG"
                    "e/a",  # The ratio between e and a of the PPG"
                    "f/a",  # The ratio between f and a of the PPG"
                    "MS",   # The PPG'(u)/PPG(peak)
                    "p2/p1",  # The PPG(p2)/PPG(p1)
                    "IPAD", # Inflection point area plus d-peak, (A2/A1) + d/a
                    "AI",   # Augmentation index is (PPG(p2) − PPG(p1))/(PPG(systolic peak) − PPG(systolic onset))
                    "SC",   # Spring constant defined as PPG"(sys) / ((SPA - MS) / SPA), derived from a physical model of the elasticity of peripheral arteries
                    "RIp1", # Reflection index using p1,(PPG(diastolic peak) − PPG(systolic onset))/(PPG(p1) − PPG(systolic onset))
                    "RIp2", # Reflection index using p2,(PPG(diastolic peak) − PPG(systolic onset))/(PPG(p2) − PPG(systolic onset))
                    ]

    df, df_features = get_features(s, fiducials, features_lst)

    return df_features

###########################################################################
######################### PPG feature extraction ##########################
###########################################################################

class features_extract_PPG:

    def __init__(self, data, peak_value, peak_time, next_peak_value, next_peak_time, onsets_values, onsets_times,
                 sample_rate,list_features, fiducials):
        """
        :param data: struct of PPG,PPG',PPG",PPG'"
            - data.sig: segment of PPG timeseries to analyse and extract features as a np array
            - data.d1: segment of PPG'
            - data.d2: segment of PPG"
            - data.d3: segment of PPG'"
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
        self.segment = data.sig
        self.segment_d1 = data.d1
        self.segment_d2 = data.d2
        self.segment_d3 = data.d3
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
        self.p1, self.p2, self.Tp1, self.Tp2 = self._getThirdDerivitivePoints()

    def map_func(self):
        """ This function assign for each name of features a function that calculates it
            :returns my_funcs: a dictionary where the key is the name of the feature and
                               the value is the function to call"""
        my_funcs = {"Tpi": self.getTpi(),
                    "Tpp": self.getTpp(),
                    "Tsys": self.getTsys(),
                    "Tdia": self.getTdia(),
                    "Tsp": self.getTsp(),
                    "Tdp": self.getTdp(),
                    "deltaT": self.get_deltaT(),
                    "Tsw10": self.getSystolicWidth_d_percent(10),
                    "Tsw25": self.getSystolicWidth_d_percent(25),
                    "Tsw33": self.getSystolicWidth_d_percent(33),
                    "Tsw50": self.getSystolicWidth_d_percent(50),
                    "Tsw66": self.getSystolicWidth_d_percent(66),
                    "Tsw75": self.getSystolicWidth_d_percent(75),
                    "Tsw90": self.getSystolicWidth_d_percent(90),
                    "Tdw10": self.getDiastolicWidth_d_percent(10),
                    "Tdw25": self.getDiastolicWidth_d_percent(25),
                    "Tdw33": self.getDiastolicWidth_d_percent(33),
                    "Tdw50": self.getDiastolicWidth_d_percent(50),
                    "Tdw66": self.getDiastolicWidth_d_percent(66),
                    "Tdw75": self.getDiastolicWidth_d_percent(75),
                    "Tdw90": self.getDiastolicWidth_d_percent(90),
                    "Tpw10": self.getSumSW_DW(10),
                    "Tpw25": self.getSumSW_DW(25),
                    "Tpw33": self.getSumSW_DW(33),
                    "Tpw50": self.getSumSW_DW(50),
                    "Tpw66": self.getSumSW_DW(66),
                    "Tpw75": self.getSumSW_DW(75),
                    "Tpw90": self.getSumSW_DW(90),
                    "Asp": self.getSystolicPeak(),
                    "Adn": self.getDicroticNotchAmplitude(),
                    "Adp": self.getDiastolicPeak(),
                    "Aoff": self.getPulseOffsetAmplitude(),
                    "AUCpi": self.getAUCpi(),
                    "AUCsys": self.getAUCsys(),
                    "AUCdia": self.getAUCdia(),
                    "IPR": self.getIPR(),
                    "Tsys/Tdia": self.get_ratio_Tsys_Tdia(),
                    "Tpw25/Tpi": self.get_ratio_Tpwx_Tpi(25),
                    "Tpw50/Tpi": self.get_ratio_Tpwx_Tpi(50),
                    "Tpw75/Tpi": self.get_ratio_Tpwx_Tpi(75),
                    "Tpw25/Tsp": self.get_ratio_Tpwx_Tsp(25),
                    "Tpw50/Tsp": self.get_ratio_Tpwx_Tsp(50),
                    "Tpw75/Tsp": self.get_ratio_Tpwx_Tsp(75),
                    "Tdw10/Tsw10": self.get_ratio_Tdwx_Tswx(10),
                    "Tdw25/Tsw25": self.get_ratio_Tdwx_Tswx(25),
                    "Tdw33/Tsw33": self.get_ratio_Tdwx_Tswx(33),
                    "Tdw50/Tsw50": self.get_ratio_Tdwx_Tswx(50),
                    "Tdw66/Tsw66": self.get_ratio_Tdwx_Tswx(66),
                    "Tdw75/Tsw75": self.get_ratio_Tdwx_Tswx(75),
                    "Tdw90/Tsw90": self.get_ratio_Tdwx_Tswx(90),
                    "Tsp/Tpi": self.get_ratio_Tsp_Tpi(),
                    "Asp/Aoff": self.get_ratio_Asp_Aoff(),
                    "Adp/Asp": self.get_ratio_Adp_Asp(),
                    "IPA": self.getIPA(),
                    "Tsp/Asp": self.get_ratio_Tsp_Asp(),
                    "Asp/deltaT": self.get_ratio_Asp_deltaT(),
                    "Asp/(Tpi-Tsp)": self.get_ratio_Asp_TpiTsp(),
                    "u": self.get_u(),
                    "v": self.get_v(),
                    "w": self.get_v(),
                    "a": self.get_a(),
                    "b": self.get_b(),
                    "c": self.get_c(),
                    "d": self.get_d(),
                    "e": self.get_e(),
                    "f": self.get_f(),
                    "Tu": self.get_Tu(),
                    "Tv": self.get_Tv(),
                    "Tw": self.get_Tw(),
                    "Ta": self.get_Ta(),
                    "Tb": self.get_Tb(),
                    "Tc": self.get_Tc(),
                    "Td": self.get_Td(),
                    "Te": self.get_Te(),
                    "Tf": self.get_Tf(),
                    "Tb–c": self.get_Tbc(),
                    "Tb–d": self.get_Tbd(),
                    "Tp1": self.get_Tp1(),
                    "Tp2": self.get_Tp2(),
                    "Tp1–dp": self.get_Tp1_dp(),
                    "Tp2–dp": self.get_Tp2_dp(),
                    "Tu/CP": self.get_ratio_Tu_CP(),
                    "Tv/CP": self.get_ratio_Tv_CP(),
                    "Tw/CP": self.get_ratio_Tw_CP(),
                    "Ta/CP": self.get_ratio_Ta_CP(),
                    "Tb/CP": self.get_ratio_Tb_CP(),
                    "Tc/CP": self.get_ratio_Tc_CP(),
                    "Td/CP": self.get_ratio_Td_CP(),
                    "Te/CP": self.get_ratio_Te_CP(),
                    "Tf/CP": self.get_ratio_Tf_CP(),
                    "(Tu-Ta)/CP": self.get_ratio_Tu_Ta_CP(),
                    "(Tv-Tb)/CP": self.get_ratio_Tv_Tb_CP(),
                    "AGI_bcdef/a": self.get_aging_index0(),
                    "AGI_bcde/a": self.get_aging_index1(),
                    "AGI_bcd/a": self.get_aging_index2(),
                    "AGI_be/a": self.get_aging_index3(),
                    "v/u": self.get_ratio_v_u(),
                    "w/u": self.get_ratio_w_u(),
                    "b/a": self.get_ratio_b_a(),
                    "c/a": self.get_ratio_c_a(),
                    "d/a": self.get_ratio_d_a(),
                    "e/a": self.get_ratio_e_a(),
                    "f/a": self.get_ratio_f_a(),
                    "IPAD": self.getIPAD(),
                    "MS": self.getMS(),
                    "AI": self.getAI(),
                    "SC": self.getSC(),
                    "p2/p1": self.get_ratio_p2_p1(),
                    "RIp1": self.getRIp1(),
                    "RIp2": self.getRIp2(),
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
        c = (self.fiducials.c-self.fiducials.os).values[0]
        d = (self.fiducials.d-self.fiducials.os).values[0]
        e = (self.fiducials.e-self.fiducials.os).values[0]
        f = (self.fiducials.f-self.fiducials.os).values[0]
        Ta = a / self.sample_rate
        Tb = b / self.sample_rate
        Tc = c / self.sample_rate
        Td = d / self.sample_rate
        Te = e / self.sample_rate
        Tf = f / self.sample_rate

        return a, b, c, d, e, f, Ta, Tb, Tc, Td, Te, Tf

    def _getThirdDerivitivePoints(self):
        """Calculate third derivitive points from a SINGLE Onset-Onset segment of PPG'"
        :return p1: Sample distance from PPG onset to the first local maximum after b on PPG'"
        :return p2: Sample distance from PPG onset to the last local minimum before d, if c = d, then the first local minimum after d on PPG'"
        :return Tp1: Time from PPG onset to the first local maximum after b on PPG'"
        :return Tp2: Time from PPG onset to last local minimum before d, if c = d, then the first local minimum after d on PPG'"
        """

        p1 = (self.fiducials.p1-self.fiducials.os).values[0]
        p2 = (self.fiducials.p2-self.fiducials.os).values[0]
        Tp1 = p1 / self.sample_rate
        Tp2 = p2 / self.sample_rate

        return p1, p2, Tp1, Tp2
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

    def getTpi(self):
        """ Tpi which means the Pulse Interval,
            the time between the pulse onset and pulse offset.
            :return Tpi feature:
        """
        return self.onsets_times[1] - self.onsets_times[0]
    def getTpp(self):
        """ Tpp means the Peak-to-Peak Interval,
            the time between two consecutive systolic peaks.
            :return  Tpp feature:
        """
        cardiac_period = self.next_peak_time - self.peak_time
        return cardiac_period

    def getTsys(self):
        """ Tsys means the Systolic Time,
            the time between the pulse onset and dicrotic notch.
        :return Tsys feature:
        """

        Tsys = self.Tdn
        return Tsys

    def getTdia(self):
        """ Tdia means the Diastolic Time,
            the time between the dicrotic notch and pulse offset.
        :return Tdia feature:
        """

        Tdia = self.getTpi()-self.Tdp
        return Tdia

    def getTsp(self):
        """ Tsp means the Systolic Peak Time,
            the time between the pulse onset and systolic peak.
        :return Tsp feature:
        """
        left_onset = self.onsets_times[0]
        Tsp = self.peak_time - left_onset
        return Tsp

    def getTdp(self):
        """ Tdp means the Diastolic Peak Time,
            the time between the pulse onset and diastolic peak.
        :return Tdp feature:
        """

        Tdp = self.Tdp
        return Tdp

    def get_deltaT(self):
        """ deltaT means the Time Delay,
            the time between the systolic peak and diastolic peak.
        :return deltaT feature:
        """

        deltaT = self.Tdp-self.getTsp()
        return deltaT

    def getSystolicWidth_d_percent(self, d):
        """ The function calculates the Systolic Width,
            the width at x% of the Systolic Peak Amplitude between the pulse onset and systolic peak.
            :param d: the percentage chosen to calculate the width
            :return Tswx feature:
        """
        # value in segment corresponding to d percent of pulse height
        d_percent_val = (d/100)*(self.peak_value - self.onsets_values[0]) + self.onsets_values[0]
        time_of_d = self._getSysTime_from_val(d_percent_val)
        Tswx = self.peak_time - (time_of_d/self.sample_rate)
        return Tswx


    def getDiastolicWidth_d_percent(self, d):
        """ The function calculates the Diastolic Width,
            the width at x% of the Systolic Peak Amplitude between the systolic peak and pulse offset.
            :param d: the percentage chosen to calculate the width
            :return Tdwx feature:
        """
        # value in segment corresponding to d percent of pulse height
        d_percent_val = (d/100)*(self.peak_value - self.onsets_values[1]) + self.onsets_values[1]
        time_of_d = self._getDiaTime_from_val(d_percent_val)
        Tdwx = (time_of_d/self.sample_rate) - self.peak_time
        return Tdwx

    def getSumSW_DW(self, d):
        """ The function calculates Pulse Width,
            the sum of the Systolic Width and Diastolic Width at x%.
            :param d: the percentage chosen to calculate the width
            :return Tpwx feature:
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        return Tswx + Tdwx

    def getSystolicPeak(self):
        """ The function calculates the Systolic Peak Amplitude,
            the difference in amplitude between the pulse onset and systolic peak.
            :return Systolic Peak Amplitude feature:
        """
        sys_peak = self.peak_value - self.onsets_values[0]
        return sys_peak
    def getDicroticNotchAmplitude(self):
        """ The function calculates the Dicrotic Notch Amplitude,
            the difference in amplitude between the pulse onset and the dicrotic notch.
            :return Dicrotic Notch Amplitude feature:
        """
        dn_value = self.segment[self.dn]
        dn_amp = dn_value - self.onsets_values[0]
        return dn_amp

    def getDiastolicPeak(self):
        """ The function calculates the Diastolic Peak Amplitude,
            the difference in amplitude between the pulse onset and the diastolic peak.
            :return Diastolic Peak Amplitude feature:
        """
        ## temp solution 04/04/2023 --> if DN==DP
        dp_value = self.segment[self.dp+20]
        dia_peak = dp_value - self.onsets_values[0]
        return dia_peak

    def getPulseOffsetAmplitude(self):
        """ The function calculates the Pulse Offset Amplitude,
            the difference in amplitude between the pulse onset and pulse offset.
            :return Pulse Offset Amplitude feature:
        """

        offset_value = self.onsets_values[1]
        Aoff = offset_value - self.onsets_values[0]
        return Aoff

    def getAUCpi(self):
        """ The function the Area Under Pulse Interval Curve,
            the area under the pulse wave between pulse onset and pulse offset.
            :return AUCpi feature:
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
        AUCpi = 10*sum/((self.peak_value - baseline) * (right_onset_time - left_onset_time))

        return AUCpi

    def getAUCsys(self):
        """ The function calculates the Area Under Systolic Curve,
            the area under the pulse wave between the pulse onset and dicrotic notch.
            :return AUCsys feature:
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
        AUCsys = 10*sum/((self.peak_value - baseline) * (right_onset_time - left_onset_time))

        return AUCsys

    def getAUCdia(self):
        """ The function calculates Area Under Diastolic Curve,
            the area under the pulse wave between the dicrotic notch and pulse offset.
            :return AUCdia feature:
        """
        AUCdia = self.getAUCpi()-self.getAUCsys()
        return AUCdia

    def getIPR(self):
        """ The function calculates the Instantaneous Pulse Rate, 60/CP.
            :return IPR feature:
        """
        IPR = 60/self.getTpp()
        return IPR

    def get_ratio_Tsys_Tdia(self):
        """ The function calculates the ratio of the Systolic Time to the Diastolic Time.
            :return Tsys/Tdia feature:
        """
        Tsys_Tdia = self.getTsys()/self.getTdia()
        return Tsys_Tdia

    def get_ratio_Tpwx_Tpi(self, d):
        """ The function calculates the ratio of the Pulse Width at x% of the Systolic Peak Amplitude to the Systolic Peak Time.
            :param d: the percentage chosen to calculate the width
            :return The ratio of the Tpi to the Pulse Width feature:
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        Tpi = self.getTpi()
        return width/Tpi

    def get_ratio_Tpwx_Tsp(self, d):
        """ The function calculates the ratio of the Pulse Width at x% of the Systolic Peak Amplitude to the Systolic Peak Time.
            :param d: the percentage chosen to calculate the width
            :return Tpwx/Tsp feature:
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        Tpwx = Tswx + Tdwx
        Tsp = self.getSystolicPeakTime()
        return Tpwx/Tsp

    def get_ratio_Tdwx_Tswx(self, d):
        """ The function calculates the ratio of the Diastolic Width to the Systolic Width at x% width.
            :param d: the percentage chosen to calculate the width
            :return Tdwx/Tswx feature:
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        return Tdwx/Tswx

    def get_ratio_Tsp_Tpi(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Pulse Interval
            :return Tsp/Tpi feature:
        """
        Tsp = self.getSystolicPeakTime()
        Tpi = self.getTpi()
        return Tsp/Tpi

    def get_ratio_Asp_Aoff(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Pulse Interval
            :return Asp/Aoff feature:
        """
        Asp = self.getSystolicPeak()
        Aoff = self.onsets_values[1]
        return Asp/Aoff

    def get_ratio_Adp_Asp(self):
        """ The function calculates Reflection Index,
            the ratio of the Diastolic Peak Amplitude to the Systolic Peak Amplitude.
            :return Reflection Index feature:
        """
        RI = self.getDiastolicPeak()/self.getSystolicPeak()
        return RI
    def getIPA(self):
        """ The function calculates the Inflection Point Area,
            the ratio of the Area Under Diastolic Curve to the Area Under Systolic Curve.
            :return IPA feature:
        """
        IPA = self.getAUCdia()/self.getAUCsys()
        return IPA

    def get_ratio_Tsp_Asp(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Systolic Peak Amplitude.
            :return Tsp/Asp feature:
        """

        return self.getTsp()/self.getSystolicPeak()
    def get_ratio_Asp_deltaT(self):
        """ The function calculates the Stiffness Index,
            the ratio of the Systolic Peak Amplitude to the Time Delay.
            :return Stiffness Index feature:
        """
        SI = self.getSystolicPeak()/(self.getTdp()-self.getTsp())
        return SI

    def get_ratio_Asp_TpiTsp(self):
        """ The function calculates the ratio of the Systolic Peak Amplitude to the difference between the Pulse Interval and Systolic Peak Time.
            :return Asp/(Tpi-Tsp):
        """
        Tpi = self.getTpi()
        Tsp = self.getSystolicPeakTime()
        Asp = self.getSystolicPeak()
        return Asp/(Tpi-Tsp)


    def get_u(self):
        """ u means the u-point sample,
            the sample between the pulse onset and u-point.
            :return u feature:
        """
        return self.u


    def get_v(self):
        """ v means the v-point sample,
            the sample between the pulse onset and v-point.
            :return v feature:
        """
        return self.v

    def get_w(self):
        """ w means the w-point sample,
            the sample between the pulse onset and w-point.
            :return w feature:
        """
        return self.w

    def get_a(self):
        """ a means the a-point sample,
            the sample between the pulse onset and a-point.
            :return a feature:
        """
        return self.a

    def get_b(self):
        """ b means the b-point sample,
            the sample between the pulse onset and b-point.
            :return b feature:
        """
        return self.b

    def get_c(self):
        """ c means the c-point sample,
            the sample between the pulse onset and c-point.
            :return c feature:
        """
        return self.c

    def get_d(self):
        """ d means the d-point sample,
            the sample between the pulse onset and d-point.
            :return d feature:
        """
        return self.d

    def get_e(self):
        """ e means the e-point sample,
            the sample between the pulse onset and e-point.
            :return e feature:
        """
        return self.e

    def get_f(self):
        """ f means the f-point sample,
            the sample between the pulse onset and f-point.
            :return f feature:
        """
        return self.f

    def get_Tu(self):
        """ Tu means the u-point time,
            the time between the pulse onset and u-point.
            :return Tu feature
        """
        return self.Tu

    def get_Tv(self):
        """ Tv means the v-point time,
            the time between the pulse onset and v-point.
            :return Tv feature:
        """
        return self.Tv

    def get_Tw(self):
        """ Tw means the v-point time,
            the time between the pulse onset and w-point.
            :return Tw feature:
        """
        return self.Tw

    def get_Ta(self):
        """ Ta means the a-point time,
            the time between the pulse onset and a-point.
            :return Ta feature:
        """
        return self.Ta

    def get_Tb(self):
        """ Tb means the b-point time,
            the time between the pulse onset and b-point.
            :return Tb feature:
        """
        return self.Tb

    def get_Tc(self):
        """ Tc means the c-point time,
            the time between the pulse onset and c-point.
            :return Tb feature:
        """
        return self.Tc

    def get_Td(self):
        """ Td means the d-point time,
            the time between the pulse onset and d-point.
            :return Td feature:
        """
        return self.Td


    def get_Te(self):
        """ Te means the e-point time,
            the time between the pulse onset and e-point.
            :return Te feature:
        """
        return self.Te

    def get_Tf(self):
        """ Tf means the f-point time,
            the time between the pulse onset and f-point.
            :return Tf feature:
        """
        return self.Tf

    def get_Tbc(self):
        """ Tbc means the b–c time,
            the time between the b-point and c-point.
            :return Tbc feature:
        """
        return self.Tc-self.Tb

    def get_Tbd(self):
        """ Tbd means the b–d time,
            the time between the b-point and d-point.
            :return Tbd feature:
        """
        return self.Td - self.Tb

    def get_Tp1(self):
        """ Tp1 means the p1-point time,
            the time between the pulse onset and p1-point.
            :return Tp1 feature:
        """
        return self.Tp1

    def get_Tp2(self):
        """ Tp2 means the p1-point time,
            the time between the pulse onset and p2-point.
            :return Tp2 feature:
        """
        return self.Tp2

    def get_Tp1_dp(self):
        """ The function calculatesthe p1–dia time,
            the time between the p1-point and diastolic peak.
            :return Tdia-Tp1 feature:
        """
        Tp1_dia=(self.dp-self.p1)/self.sample_rate
        return Tp1_dia

    def get_Tp2_dp(self):
        """ The function calculatesthe p2–dia time,
            the time between the p2-point and diastolic peak.
            :return Tdia-Tp2 feature:
        """
        Tp2_dia=(self.dp-self.p2)/self.sample_rate
        return Tp2_dia

    
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

    def getRIp1(self):
        """ The function calculates Reflection index using p1, (PPG(diastolic peak) − PPG(systolic onset))/(PPG(p1) − PPG(systolic onset))
            :return RIp1 feature:
        """
        RIp1 = self.getDiastolicPeak()/self.segment[self.p1]
        return RIp1

    def getRIp2(self):
        """ The function calculates Reflection index using p2, (PPG(diastolic peak) − PPG(systolic onset))/(PPG(p2) − PPG(systolic onset))
            :return RIp2 feature:
        """
        RIp2 = self.getDiastolicPeak()/self.segment[self.p2]
        return RIp2

    def getAI(self):
        """ The function calculates the Augmentation index, which is (x(p2) − x(p1))/(x(sys) − x(0)), where is x the amplitude
            :return pir feature:
        """
        AI = (self.segment[self.p2]-self.segment[self.p1])/self.peak_value
        return AI

    def getIPAD(self):
        """ The function calculates the inflection point area plus d-peak, (A2/A1)+PPG"(d)/PPG"(a)
            :return IPAD feature:
        """
        IPAD = self.getAUCdia()/self.getAUCsys()+self.get_ratio_d_a()
        return IPAD

    def getMS(self):
        """ The function calculates Maximum slope, PPG'(u)/(PPG(systolic peak) − PPG(systolic onset))
            :return MS feature:
        """
        MS = self.segment_d1[self.u]
        return MS

    def getSC(self):
        """ The function calculates Spring constant, PPG"(systolic peak)/((SPA-MS)/SPA),
            :return SC feature:
        """
        ddxSPA=self.segment_d2[(self.getTsp()*self.sample_rate).astype(int)]
        SPA=self.getSystolicPeak()
        MS=self.segment[self.u]
        SC = ddxSPA/((SPA-MS)/SPA)
        return SC

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

    def get_ratio_Tu_CP(self):
        """ The function calculates the ratio of Tu and CP.
            :return Tu and CP ratio:
        """
        T1 = self.get_Tu()
        return T1 / self.getTpp()


    def get_ratio_Tv_CP(self):
        """ The function calculates the ratio of Tv and CP.
            :return Tv and CP ratio:
        """
        Tv = self.get_Tv()
        return Tv / self.getTpp()

    def get_ratio_Tw_CP(self):
        """ The function calculates the ratio of Tw and CP.
            :return Tw and CP ratio:
        """
        Tw = self.get_Tw()
        return Tw / self.getTpp()

    def get_ratio_Ta_CP(self):
        """ The function calculates the ratio of Ta and CP.
            :return Ta and CP ratio feature:
        """
        Ta = self.get_Ta()
        return Ta / self.getTpp()

    def get_ratio_Tb_CP(self):
        """ The function calculates the ratio of Tb and CP.
            :return Tb and CP ratio feature:
        """
        Tb = self.get_Tb()
        return Tb / self.getTpp()

    def get_ratio_Tc_CP(self):
        """ The function calculates the ratio of Tc and CP.
            :return Tc and CP ratio feature:
        """
        Tc = self.get_Tc()
        return Tc / self.getTpp()

    def get_ratio_Td_CP(self):
        """ The function calculates the ratio of Td and CP.
            :return Td and CP ratio feature:
        """
        Td = self.get_Td()
        return Td / self.getTpp()

    def get_ratio_Te_CP(self):
        """ The function calculates the ratio of Te and CP.
            :return Te and CP ratio feature:
        """
        Te = self.get_Te()
        return Te / self.getTpp()

    def get_ratio_Tf_CP(self):
        """ The function calculates the ratio of Tf and CP.
            :return Tf and CP ratio feature:
        """
        Tf = self.get_Tf()
        return Tf / self.getTpp()

    def get_ratio_p2_p1(self):
        """ The function calculates the PPG(p2)/PPG(p1).
            :return p2 and p1 ratio feature:
        """
        Rp2p1 = self.segment[self.p2]/self.segment[self.p1]
        return Rp2p1

    def get_ratio_Tu_Ta_CP(self):
        """ The function calculates the ratio between the interval maximum/minimum peaks of 1st derivative and CP
            :return (Tu - Ta) / CP:
        """
        Ta = self.get_Ta()
        Tu = self.get_Tu()
        CP = self.getTpp()
        return (Tu - Ta) / CP

    def get_ratio_Tv_Tb_CP(self):
        """ The function calculates the ratio between the interval
            :return (Tv - Tb) / CP:
        """
        Tb = self.get_Tb()
        return (self.get_Tv() - Tb) / self.getTpp()

    def get_aging_index0(self):
        """ The function calculates the aging index of (PPG"(b)-PPG"(c)-PPG"(d)-PPG"(e)-PPG"(f))/PPG"(a).
            :return (B-C-D-E-F)/A:
        """
        A = self.segment_d2[self.get_a()]
        B = self.segment_d2[self.get_b()]
        C = self.segment_d2[self.get_c()]
        D = self.segment_d2[self.get_d()]
        E = self.segment_d2[self.get_e()]
        F = self.segment_d2[self.get_f()]
        return (B-C-D-E-F)/A
    def get_aging_index1(self):
        """ The function calculates the aging index of (PPG"(b)-PPG"(c)-PPG"(d)-PPG"(e))/PPG"(a).
            :return (B-C-D-E)/A:
        """
        A = self.segment_d2[self.get_a()]
        B = self.segment_d2[self.get_b()]
        C = self.segment_d2[self.get_c()]
        D = self.segment_d2[self.get_d()]
        E = self.segment_d2[self.get_e()]
        return (B-C-D-E)/A

    def get_aging_index2(self):
        """ The function calculates the aging index of (PPG"(b)-PPG"(c)-PPG"(d))/PPG"(a).
            :return (B-C-D)/A:
        """
        A = self.segment_d2[self.get_a()]
        B = self.segment_d2[self.get_b()]
        C = self.segment_d2[self.get_c()]
        D = self.segment_d2[self.get_d()]
        return (B-C-D)/A

    def get_aging_index3(self):
        """ The function calculates the aging index of (PPG"(b)-PPG"(e))/PPG"(a).
            :return (B-E)/A:
        """
        A = self.segment_d2[self.get_a()]
        B = self.segment_d2[self.get_b()]
        E = self.segment_d2[self.get_e()]
        return (B-E)/A

    def get_ratio_v_u(self):
        """ This function calculates PPG'(v)/PPG'(u).
            :return u_max and u_max ratio:
        """
        u_max = self.segment_d1[self.get_u()]
        v_min = self.segment_d1[self.get_v()]
        return v_min / u_max

    def get_ratio_w_u(self):
        """ This function calculates PPG'(w)/PPG'(u).
            :return u_max and u_max ratio:
        """
        u_max = self.segment_d1[self.get_u()]
        w_max = self.segment_d1[self.get_w()]
        return w_max / u_max

    def get_ratio_b_a(self):
        """ This function calculates PPG"(b)/PPG"(a).
            :return b_min and a_max ratio feature:
        """
        a_max = self.segment_d2[self.get_a()]
        b_min = self.segment_d2[self.get_b()]
        return b_min / a_max

    def get_ratio_c_a(self):
        """ This function calculates PPG"(c)/PPG"(a).
            :return c_max and a_max ratio feature:
        """
        a_max = self.segment_d2[self.get_a()]
        c_max = self.segment_d2[self.get_c()]
        return c_max / a_max

    def get_ratio_d_a(self):
        """ This function calculates PPG"(d)/PPG"(a).
            :return d_min and a_max ratio feature:
        """
        a_max = self.segment_d2[self.get_a()]
        d_min = self.segment_d2[self.get_d()]
        return d_min / a_max

    def get_ratio_e_a(self):
        """ This function calculates PPG"(e)/PPG"(a).
            :return e_max and a_max ratio feature:
        """
        a_max = self.segment_d2[self.get_a()]
        e_max = self.segment_d2[self.get_e()]
        return e_max / a_max

    def get_ratio_f_a(self):
        """ This function calculates PPG"(f)/PPG"(a).
            :return f_min and a_max ratio feature:
        """
        a_max = self.segment_d2[self.get_a()]
        f_min = self.segment_d2[self.get_f()]
        return f_min / a_max

###########################################################################
############################# Get PPG features ############################
###########################################################################

def get_features(s, fiducials, features_lst):
    """
    The function calculates the biomedical features of PPG signal.

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param fiducials: M-d Dateframe, where M is the number of fiducial points
    :param features_lst: list of features

    :return
        - df: data frame with onsets, offset and peaks
        - df_features: data frame with PPG signal features
    """

    fs=s.fs
    ppg=s.filt_sig
    data = DotMap()

    df = pd.DataFrame()
    df_features = pd.DataFrame(columns=features_lst)
    peaks = fiducials.pk.values
    onsets = fiducials.os.values

    # display(df_features)
    for i in range(len(onsets) - 1):
        #         #     print(f'i is {i}')
        onset = onsets[i]
        offset = onsets[i + 1]
        data.sig = ppg[int(onset):int(offset)]
        data.d1 = s.filt_d1[int(onset):int(offset)]
        data.d2 = s.filt_d2[int(onset):int(offset)]
        data.d3 = s.filt_d3[int(onset):int(offset)]
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
            #         plt.plot(data.sig)
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