import pyPPG

import numpy as np
from dotmap import DotMap
import pandas as pd
from scipy.signal import find_peaks

###########################################################################
####################### PPG biomarkers extraction #########################
###########################################################################

class BmExctator:
    """
    Class that extracts the PPG biomarkers.
    """

    def __init__(self, data: DotMap, peak_value: float, peak_time: float, next_peak_value: float, next_peak_time: float, onsets_values: np.array, onsets_times: np.array,
                 sample_rate: int, list_biomarkers: list, fiducials: pd.DataFrame):
        """

        :param data: struct of PPG,PPG', PPG", PPG'":

            - data.sig: segment of PPG timeseries to analyse and extract biomarkers as a np array
            - data.d1: segment of PPG'
            - data.d2: segment of PPG"
            - data.d3: segment of PPG'"
        :type data: DotMap
        :param peak_value: PPG peak value
        :type peak_value: float
        :param peak_time: the time corresponding to the peak detected
        :type peak_time: float
        :param next_peak_value: PPG next peak value
        :type next_peak_value: float
        :param next_peak_time: the time corresponding to the peak detected
        :type next_peak_time: float
        :param onsets_values: array of PPG two onsets values surrounding the peak
        :type onsets_values: Numpy array
        :param onsets_times: array of the two times corresponding to each onset detected
        :type onsets_times: Numpy array
        :param sample_rate: segment data sample rate
        :type sample_rate: int
        :param list_biomarkers: list of biomarkers
        :type list_biomarkers: list
        :param fiducials: location of fiducial points of the given pulse wave
        :type fiducials: DataFrame

        """

        self.fiducials=fiducials
        self.list_biomarkers=list_biomarkers
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
        """ This function assign for each name of biomarkers a function that calculates it

        :returns my_funcs: a dictionary where the key is the name of the biomarker and the value is the function to call
        """
        my_funcs = {#PPG signal
                    "Tpi": self.getTpi(),
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

                    # Signal ratios
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

                    # PPG derivatives
                    "u": self.get_u(),
                    "v": self.get_v(),
                    "w": self.get_w(),
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
                    "Tb-c": self.get_Tbc(),
                    "Tb-d": self.get_Tbd(),
                    "Tp1": self.get_Tp1(),
                    "Tp2": self.get_Tp2(),
                    "Tp1-dp": self.get_Tp1_dp(),
                    "Tp2-dp": self.get_Tp2_dp(),

                    # Derivatives ratios
                    "Tu/Tpi": self.get_ratio_Tu_Tpi(),
                    "Tv/Tpi": self.get_ratio_Tv_Tpi(),
                    "Tw/Tpi": self.get_ratio_Tw_Tpi(),
                    "Ta/Tpi": self.get_ratio_Ta_Tpi(),
                    "Tb/Tpi": self.get_ratio_Tb_Tpi(),
                    "Tc/Tpi": self.get_ratio_Tc_Tpi(),
                    "Td/Tpi": self.get_ratio_Td_Tpi(),
                    "Te/Tpi": self.get_ratio_Te_Tpi(),
                    "Tf/Tpi": self.get_ratio_Tf_Tpi(),
                    "(Tu-Ta)/Tpi": self.get_ratio_TuTa_Tpi(),
                    "(Tv-Tb)/Tpi": self.get_ratio_TvTb_Tpi(),
                    "Au/Asp":self.get_ratio_Au_Asp(),
                    "Av/Au":self.get_ratio_Av_Au(),
                    "Aw/Au":self.get_ratio_Aw_Au(),
                    "Ab/Aa":self.get_ratio_Ab_Aa(),
                    "Ac/Aa":self.get_ratio_Ac_Aa(),
                    "Ad/Aa":self.get_ratio_Ad_Aa(),
                    "Ae/Aa":self.get_ratio_Ae_Aa(),
                    "Af/Aa":self.get_ratio_Af_Aa(),
                    "Ap2/Ap1": self.get_ratio_Ap2_Ap1(),
                    "(Ac-Ab)/Aa":self.get_ratio_AcAb_Aa(),
                    "(Ad-Ab)/Aa":self.get_ratio_AdAb_Aa(),
                    "AGI":self.getAGI(),
                    "AGImod":self.getAGImod(),
                    "AGIinf":self.getAGIinf(),
                    "AI": self.getAI(),
                    "RIp1": self.getRIp2(),
                    "RIp2": self.getRIp2(),
                    "SC": self.getSC(),
                    "IPAD": self.getIPAD(),
        }
        return my_funcs

    def get_biomarker_extract_func(self):
        """ This function go through the list of biomarkers and call the function that is relevant to calculate it each biomarker takes two spots in the biomarker vector: one for the average value and one for its standard deviation.

            :returns biomarkers_vec: a vector of each biomarker avg value and std of the patient"""
        my_funcs = self.map_func()
        biomarkers_vec = []

        for biomarker in self.list_biomarkers:
            func_to_call = my_funcs[biomarker]
            biomarkers_vec.append(func_to_call)
        return biomarkers_vec

    def _getPeaksOnsets(self,x):
        """Find the peaks and onsets of a short FILTERED segment of PPG

        :return: peaks
        :return: onsets
        """
        peaks, _ = find_peaks(x)
        onsets, _ = find_peaks(-x)
        return peaks, onsets

    def _getDicroticNotchDiastolicPeak(self):
        """Calculate Dicrotic Notch and Diastolic Peak of PPG

        :return: dn: Sample distance from PPG onset to the Dicrotic Notch on PPG
        :return: dp: Sample distance from PPG onset to the Diastolic Peak on PPG
        :return: Tdn: Time from PPG onset to the Dicrotic Notch on PPG
        :return: Tdp: Time from PPG onset to the Diastolic Peak on PPG
        """

        dn = (self.fiducials.dn-self.fiducials.on).values[0]
        dp = (self.fiducials.dp-self.fiducials.on).values[0]
        Tdn = dn / self.sample_rate
        Tdp = dp / self.sample_rate

        return dn, dp, Tdn, Tdp
    def _getFirstDerivitivePoints(self):
        """Calculate first derivitive points from a SINGLE Onset-Onset segment of PPG'

        :return: u: Sample distance from PPG onset to the greatest maximum peak between the left systolic onset and the systolic peak on PPG'
        :return: v: Sample distance from PPG onset to the lowest minimum pits between the systolic peak and the right systolic onset on PPG'
        :return: w: Sample distance from PPG onset to the first maximum peak after v on PPG'
        :return: Tu: Time from PPG onset to the greatest maximum peak between the left systolic onset and the systolic peak on PPG'
        :return: Tv: Time from PPG onset to the lowest minimum pits between the systolic peak and the right systolic onset on PPG'
        :return: Tw: Time from PPG onset to the first maximum peak after v on PPG'
        """

        u = (self.fiducials.u-self.fiducials.on).values[0]
        v = (self.fiducials.v-self.fiducials.on).values[0]
        w = (self.fiducials.w - self.fiducials.on).values[0]
        Tu = u / self.sample_rate
        Tv = v / self.sample_rate
        Tw = w / self.sample_rate

        return u, v, w, Tu, Tv, Tw

    def _getSecondDerivitivePoints(self):
        """Calculate second derivitive points from a SINGLE Onset-Onset segment of PPG"

        :return: a: Sample distance from PPG onset to the first maximum peak between left systolic onset and systolic peak on PPG"
        :return: b: Sample distance from PPG onset to the first minimum pits after a on PPG"
        :return: c: Sample distance from PPG onset to the greatest maximum peak between b and e, or if no maximum peak is present then the inflection point on PPG"
        :return: d: Sample distance from PPG onset to the lowest minimum pits between c and e, or if no minimum pits is present then the inflection point on PPG"
        :return: e: Sample distance from PPG onset to the greatest maximum peak between the systolic peak and  the right systolic onset on PPG"
        :return: f: Sample distance from PPG onset to the first minimum pits after e on PPG"
        :return: Ta: Time from PPG onset to the first maximum peak between left systolic onset and systolic peak on PPG"
        :return: Tb: Time from PPG onset to the first minimum pits after a on PPG"
        :return: Tc: Time from PPG onset to the greatest maximum peak between b and e, or if no maximum peak is present then the inflection point on PPG"
        :return: Td: Time from PPG onset to the lowest minimum pits between c and e, or if no minimum pits is present then the inflection point on PPG"
        :return: Te: Time from PPG onset to the greatest maximum peak between the systolic peak and  the right systolic onset on PPG"
        :return: Tf: Time from PPG onset to the first minimum pits after e on PPG"
        """

        a = (self.fiducials.a-self.fiducials.on).values[0]
        b = (self.fiducials.b-self.fiducials.on).values[0]
        c = (self.fiducials.c-self.fiducials.on).values[0]
        d = (self.fiducials.d-self.fiducials.on).values[0]
        e = (self.fiducials.e-self.fiducials.on).values[0]
        f = (self.fiducials.f-self.fiducials.on).values[0]
        Ta = a / self.sample_rate
        Tb = b / self.sample_rate
        Tc = c / self.sample_rate
        Td = d / self.sample_rate
        Te = e / self.sample_rate
        Tf = f / self.sample_rate

        return a, b, c, d, e, f, Ta, Tb, Tc, Td, Te, Tf

    def _getThirdDerivitivePoints(self):
        """Calculate third derivitive points from a SINGLE Onset-Onset segment of PPG'"

        :return: p1: Sample distance from PPG onset to the first local maximum after b on PPG'"
        :return: p2: Sample distance from PPG onset to the last local minimum before d, if c = d, then the first local minimum after d on PPG'"
        :return: Tp1: Time from PPG onset to the first local maximum after b on PPG'"
        :return: Tp2: Time from PPG onset to last local minimum before d, if c = d, then the first local minimum after d on PPG'"
        """

        p1 = (self.fiducials.p1-self.fiducials.on).values[0]
        p2 = (self.fiducials.p2-self.fiducials.on).values[0]
        Tp1 = p1 / self.sample_rate
        Tp2 = p2 / self.sample_rate

        return p1, p2, Tp1, Tp2
    def _find_nearest(self, arr, value):
        """ This function calculates the index in an array of the closest value to the arg value

            :param arr: the array where to find the index
            :param value: the value to be compared with

            :return: idx: the index of the value closest to the arg value
            """
        idx = (np.abs(arr - value)).argmin()
        return idx

    def _getTime(self, vec, val):
        """ get the time of a value in the PPG waveform  data vector
            :param vec: the array where to find the index
            :param val: the value to be compared with
            :return: index: the index of the value closest to the arg value
            """
        tmp_vec = np.array([vec[i]-val for i in range(0, len(vec))])
        index = self._find_nearest(tmp_vec, 0)
        return index

    def _getSysTime_from_val(self, val):
        """ get the time of a value in the PPG waveform  data vector

            :param val: the value from which we need the time in the timeserie

            :return: t_data: the time corresponding to the value
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

            :return: t_data: the time corresponding to the value
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

            :return: baseline slope
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

            :return: cst
        """
        left_onset_value = self.onsets_values[0]
        right_onset_value = self.onsets_values[1]
        cst = right_onset_value - left_onset_value
        return cst

    def getTpi(self):
        """ Tpi which means the Pulse Interval, the time between the pulse onset and pulse offset.

            :return: Tpi biomarker
        """
        return self.onsets_times[1] - self.onsets_times[0]
    def getTpp(self):
        """ Tpp means the Peak-to-Peak Interval, the time between two consecutive systolic peaks.

            :return:  Tpp biomarker
        """
        cardiac_period = self.next_peak_time - self.peak_time
        return cardiac_period

    def getTsys(self):
        """ Tsys means the Systolic Time, the time between the pulse onset and dicrotic notch.

        :return: Tsys biomarker
        """

        Tsys = self.Tdn
        return Tsys

    def getTdia(self):
        """ Tdia means the Diastolic Time, the time between the dicrotic notch and pulse offset.

        :return: Tdia biomarker
        """

        Tdia = self.getTpi()-self.Tdp
        return Tdia

    def getTsp(self):
        """ Tsp means the Systolic Peak Time, the time between the pulse onset and systolic peak.

        :return: Tsp biomarker
        """
        left_onset = self.onsets_times[0]
        Tsp = self.peak_time - left_onset
        return Tsp

    def getTdp(self):
        """ Tdp means the Diastolic Peak Time, the time between the pulse onset and diastolic peak.

        :return: Tdp biomarker
        """

        Tdp = self.Tdp
        return Tdp

    def get_deltaT(self):
        """ deltaT means the Time Delay, the time between the systolic peak and diastolic peak.

        :return: deltaT biomarker
        """

        deltaT = self.Tdp-self.getTsp()
        return deltaT

    def getSystolicWidth_d_percent(self, d):
        """ The function calculates the Systolic Width,
            the width at x% of the Systolic Peak Amplitude between the pulse onset and systolic peak.

            :param d: the percentage chosen to calculate the width

            :return: Tswx biomarker
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

            :return: Tdwx biomarker
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

            :return: Tpwx biomarker
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        return Tswx + Tdwx

    def getSystolicPeak(self):
        """ The function calculates the Systolic Peak Amplitude, the difference in amplitude between the pulse onset and systolic peak.

            :return: Systolic Peak Amplitude biomarker
        """
        sys_peak = self.peak_value - self.onsets_values[0]
        return sys_peak
    def getDicroticNotchAmplitude(self):
        """ The function calculates the Dicrotic Notch Amplitude, the difference in amplitude between the pulse onset and the dicrotic notch.

            :return: Dicrotic Notch Amplitude biomarker
        """
        dn_value = self.segment[self.dn]
        dn_amp = dn_value - self.onsets_values[0]
        return dn_amp

    def getDiastolicPeak(self):
        """ The function calculates the Diastolic Peak Amplitude, the difference in amplitude between the pulse onset and the diastolic peak.

            :return: Diastolic Peak Amplitude biomarker
        """
        ## temp solution 04/04/2023 --> if DN==DP
        dp_value = self.segment[self.dp+20]
        dia_peak = dp_value - self.onsets_values[0]
        return dia_peak

    def getPulseOffsetAmplitude(self):
        """ The function calculates the Pulse Offset Amplitude, the difference in amplitude between the pulse onset and pulse offset.

            :return: Pulse Offset Amplitude biomarker
        """

        offset_value = self.onsets_values[1]
        Aoff = offset_value - self.onsets_values[0]
        return Aoff

    def getAUCpi(self):
        """ The function the Area Under Pulse Interval Curve, the area under the pulse wave between pulse onset and pulse offset.

            :return: AUCpi biomarker
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
        """ The function calculates the Area Under Systolic Curve, the area under the pulse wave between the pulse onset and dicrotic notch.

            :return: AUCsys biomarker
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
        """ The function calculates Area Under Diastolic Curve, the area under the pulse wave between the dicrotic notch and pulse offset.

            :return: AUCdia biomarker
        """
        AUCdia = self.getAUCpi()-self.getAUCsys()
        return AUCdia

    def getIPR(self):
        """ The function calculates the Instantaneous Pulse Rate, 60/CP.

            :return: IPR biomarker
        """
        IPR = 60/self.getTpp()
        return IPR

    def get_ratio_Tsys_Tdia(self):
        """ The function calculates the ratio of the Systolic Time to the Diastolic Time.

            :return: Tsys/Tdia biomarker
        """
        Tsys_Tdia = self.getTsys()/self.getTdia()
        return Tsys_Tdia

    def get_ratio_Tpwx_Tpi(self, d):
        """ The function calculates the ratio of the Pulse Width at x% of the Systolic Peak Amplitude to the Systolic Peak Time.

            :param d: the percentage chosen to calculate the width

            :return: The ratio of the Tpi to the Pulse Width biomarker
        """
        sys_width = self.getSystolicWidth_d_percent(d)
        dia_width = self.getDiastolicWidth_d_percent(d)
        width = sys_width + dia_width
        Tpi = self.getTpi()
        return width/Tpi

    def get_ratio_Tpwx_Tsp(self, d):
        """ The function calculates the ratio of the Pulse Width at x% of the Systolic Peak Amplitude to the Systolic Peak Time.

            :param d: the percentage chosen to calculate the width

            :return: Tpwx/Tsp biomarker
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        Tpwx = Tswx + Tdwx
        Tsp = self.getSystolicPeakTime()
        return Tpwx/Tsp

    def get_ratio_Tdwx_Tswx(self, d):
        """ The function calculates the ratio of the Diastolic Width to the Systolic Width at x% width.

            :param d: the percentage chosen to calculate the width

            :return: Tdwx/Tswx biomarker
        """
        Tswx = self.getSystolicWidth_d_percent(d)
        Tdwx = self.getDiastolicWidth_d_percent(d)
        return Tdwx/Tswx

    def get_ratio_Tsp_Tpi(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Pulse Interval

            :return: Tsp/Tpi biomarker
        """
        Tsp = self.getSystolicPeakTime()
        Tpi = self.getTpi()
        return Tsp/Tpi

    def get_ratio_Asp_Aoff(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Pulse Interval

            :return: Asp/Aoff biomarker
        """
        Asp = self.getSystolicPeak()
        Aoff = self.onsets_values[1]
        return Asp/Aoff

    def get_ratio_Adp_Asp(self):
        """ The function calculates Reflection Index, the ratio of the Diastolic Peak Amplitude to the Systolic Peak Amplitude.

            :return: Reflection Index biomarker
        """
        RI = self.getDiastolicPeak()/self.getSystolicPeak()
        return RI
    def getIPA(self):
        """ The function calculates the Inflection Point Area, the ratio of the Area Under Diastolic Curve to the Area Under Systolic Curve.

            :return: IPA biomarker
        """
        IPA = self.getAUCdia()/self.getAUCsys()
        return IPA

    def get_ratio_Tsp_Asp(self):
        """ The function calculates the ratio of the Systolic Peak Time to the Systolic Peak Amplitude.

            :return: Tsp/Asp biomarker
        """

        return self.getTsp()/self.getSystolicPeak()
    def get_ratio_Asp_deltaT(self):
        """ The function calculates the Stiffness Index, the ratio of the Systolic Peak Amplitude to the Time Delay.

            :return: Stiffness Index biomarker
        """
        SI = self.getSystolicPeak()/(self.getTdp()-self.getTsp())
        return SI

    def get_ratio_Asp_TpiTsp(self):
        """ The function calculates the ratio of the Systolic Peak Amplitude to the difference between the Pulse Interval and Systolic Peak Time.

            :return: Asp/(Tpi-Tsp)
        """
        Tpi = self.getTpi()
        Tsp = self.getSystolicPeakTime()
        Asp = self.getSystolicPeak()
        return Asp/(Tpi-Tsp)


    def get_u(self):
        """ u means the u-point sample, the sample between the pulse onset and u-point.

            :return: u biomarker
        """
        return self.u


    def get_v(self):
        """ v means the v-point sample, the sample between the pulse onset and v-point.

            :return: v biomarker
        """
        return self.v

    def get_w(self):
        """ w means the w-point sample, the sample between the pulse onset and w-point.

            :return: w biomarker
        """
        return self.w

    def get_a(self):
        """ a means the a-point sample, the sample between the pulse onset and a-point.

            :return: a biomarker
        """
        return self.a

    def get_b(self):
        """ b means the b-point sample, the sample between the pulse onset and b-point.

            :return: b biomarker
        """
        return self.b

    def get_c(self):
        """ c means the c-point sample, the sample between the pulse onset and c-point.

            :return: c biomarker
        """
        return self.c

    def get_d(self):
        """ d means the d-point sample, the sample between the pulse onset and d-point.

            :return: d biomarker
        """
        return self.d

    def get_e(self):
        """ e means the e-point sample, the sample between the pulse onset and e-point.

            :return: e biomarker
        """
        return self.e

    def get_f(self):
        """ f means the f-point sample, the sample between the pulse onset and f-point.

            :return: f biomarker
        """
        return self.f

    def get_Tu(self):
        """ Tu means the u-point time, the time between the pulse onset and u-point.

            :return: Tu biomarker
        """
        return self.Tu

    def get_Tv(self):
        """ Tv means the v-point time, the time between the pulse onset and v-point.

            :return: Tv biomarker
        """
        return self.Tv

    def get_Tw(self):
        """ Tw means the v-point time, the time between the pulse onset and w-point.

            :return: Tw biomarker
        """
        return self.Tw

    def get_Ta(self):
        """ Ta means the a-point time, the time between the pulse onset and a-point.

            :return: Ta biomarker
        """
        return self.Ta

    def get_Tb(self):
        """ Tb means the b-point time, the time between the pulse onset and b-point.

            :return: Tb biomarker
        """
        return self.Tb

    def get_Tc(self):
        """ Tc means the c-point time, the time between the pulse onset and c-point.

            :return: Tb biomarker
        """
        return self.Tc

    def get_Td(self):
        """ Td means the d-point time, the time between the pulse onset and d-point.

            :return: Td biomarker
        """
        return self.Td


    def get_Te(self):
        """ Te means the e-point time, the time between the pulse onset and e-point.

            :return: Te biomarker
        """
        return self.Te

    def get_Tf(self):
        """ Tf means the f-point time, the time between the pulse onset and f-point.

            :return: Tf biomarker
        """
        return self.Tf

    def get_Tbc(self):
        """ Tbc means the b–c time, the time between the b-point and c-point.

            :return: Tbc biomarker
        """
        return self.Tc-self.Tb

    def get_Tbd(self):
        """ Tbd means the b–d time, the time between the b-point and d-point.

            :return: Tbd biomarker
        """
        return self.Td - self.Tb

    def get_Tp1(self):
        """ Tp1 means the p1-point time, the time between the pulse onset and p1-point.

            :return: Tp1 biomarker
        """
        return self.Tp1

    def get_Tp2(self):
        """ Tp2 means the p1-point time, the time between the pulse onset and p2-point.

            :return: Tp2 biomarker
        """
        return self.Tp2

    def get_Tp1_dp(self):
        """ The function calculatesthe p1–dia time, the time between the p1-point and diastolic peak.

            :return: Tdia-Tp1 biomarker
        """
        Tp1_dia=(self.dp-self.p1)/self.sample_rate
        return Tp1_dia

    def get_Tp2_dp(self):
        """ The function calculatesthe p2–dia time, the time between the p2-point and diastolic peak.

            :return: Tdia-Tp2 biomarker
        """
        Tp2_dia=(self.dp-self.p2)/self.sample_rate
        return Tp2_dia


    def get_ratio_Tu_Tpi(self):
        """ The function calculates the ratio of the u-point time to the Pulse Interval.

            :return: Tu/Tpi biomarker
        """
        T1 = self.get_Tu()
        return T1 / self.getTpp()


    def get_ratio_Tv_Tpi(self):
        """ The function calculates the ratio of the v-point time to the Pulse Interval.

            :return: Tv/Tpi biomarker
        """
        Tv = self.get_Tv()
        return Tv / self.getTpp()

    def get_ratio_Tw_Tpi(self):
        """ The function calculates the ratio of the w-point time to the Pulse Interval.

            :return: Tw/Tpi biomarker
        """
        Tw = self.get_Tw()
        return Tw / self.getTpp()

    def get_ratio_Ta_Tpi(self):
        """ The function calculates the ratio of the a-point time to the Pulse Interval.

            :return: Ta/Tpi biomarker
        """
        Ta = self.get_Ta()
        return Ta / self.getTpp()

    def get_ratio_Tb_Tpi(self):
        """ The function calculates the ratio of the b-point time to the Pulse Interval.

            :return: Tb/Tpi biomarker
        """
        Tb = self.get_Tb()
        return Tb / self.getTpp()

    def get_ratio_Tc_Tpi(self):
        """ The function calculates the ratio of the c-point time to the Pulse Interval.

            :return: Tc/Tpi biomarker
        """
        Tc = self.get_Tc()
        return Tc / self.getTpp()

    def get_ratio_Td_Tpi(self):
        """ The function calculates the ratio of the d-point time to the Pulse Interval.

            :return: Td/Tpi biomarker
        """
        Td = self.get_Td()
        return Td / self.getTpp()

    def get_ratio_Te_Tpi(self):
        """ The function calculates the ratio of the e-point time to the Pulse Interval.

            :return: Te/Tpi biomarker
        """
        Te = self.get_Te()
        return Te / self.getTpp()

    def get_ratio_Tf_Tpi(self):
        """ The function calculates the ratio of the f-point time to the Pulse Interval.

            :return: Tf/Tpi biomarker
        """
        Tf = self.get_Tf()
        return Tf / self.getTpp()

    def get_ratio_TuTa_Tpi(self):
        """ The function calculates the ratio of the difference between the u-point time and a-point time to the Pulse Interval.

            :return: (Tu-Ta)/Tpi biomarker
        """
        Ta = self.get_Ta()
        Tu = self.get_Tu()
        Tpi = self.getTpi()
        return (Tu - Ta) / Tpi

    def get_ratio_TvTb_Tpi(self):
        """ The function calculates the ratio of the difference between the v-point time and b-point time to the Pulse Interval.

            :return: (Tv-Tb)/Tpi
        """
        Tb = self.get_Tb()
        Tv = self.get_Tv()
        Tpi = self.getTpi()
        return (Tv-Tb)/Tpi

    def get_ratio_Au_Asp(self):
        """ This function calculates the ratio of the u-point amplitude to the Systolic Peak Amplitude.

            :return: Au/Asp biomarker
        """
        u_max = self.segment_d1[self.get_u()]
        sp_amp = self.peak_value
        return u_max / sp_amp

    def get_ratio_Av_Au(self):
        """ This function calculates the ratio of the v-point amplitude to the u-point amplitude.

            :return: Av/Au biomarker
        """
        u_max = self.segment_d1[self.get_u()]
        v_min = self.segment_d1[self.get_v()]
        return v_min / u_max

    def get_ratio_Aw_Au(self):
        """ This function calculates the ratio of the w-point amplitude to the u-point amplitude.

            :return: Aw/Au biomarker
        """
        u_max = self.segment_d1[self.get_u()]
        w_max = self.segment_d1[self.get_w()]
        return w_max / u_max

    def get_ratio_Ab_Aa(self):
        """ This function calculates the ratio of the b-point amplitude to the a-point amplitude.

            :return: Ab/Aa biomarker
        """
        a_max = self.segment_d2[self.get_a()]
        b_min = self.segment_d2[self.get_b()]
        return b_min / a_max

    def get_ratio_Ac_Aa(self):
        """ This function calculates the ratio of the c-point amplitude to the a-point amplitude.

            :return: Ac/Aa biomarker
        """
        a_max = self.segment_d2[self.get_a()]
        c_max = self.segment_d2[self.get_c()]
        return c_max / a_max

    def get_ratio_Ad_Aa(self):
        """ This function calculates the ratio of the d-point amplitude to the a-point amplitude.

            :return: Ad/Aa biomarker
        """
        a_max = self.segment_d2[self.get_a()]
        d_min = self.segment_d2[self.get_d()]
        return d_min / a_max

    def get_ratio_Ae_Aa(self):
        """ This function calculates the ratio of the e-point amplitude to the a-point amplitude.

            :return: Ae/Aa biomarker
        """
        a_max = self.segment_d2[self.get_a()]
        e_max = self.segment_d2[self.get_e()]
        return e_max / a_max

    def get_ratio_Af_Aa(self):
        """ This function calculates the ratio of the f-point amplitude to the a-point amplitude.

            :return: Af/Aa biomarker
        """
        a_max = self.segment_d2[self.get_a()]
        f_min = self.segment_d2[self.get_f()]
        return f_min / a_max

    def get_ratio_Ap2_Ap1(self):
        """ The function calculates the ratio of the p2-point amplitude to the p1-point amplitude.

            :return: Ap2/Ap1 biomarker
        """
        Rp2p1 = self.segment[self.p2]/self.segment[self.p1]
        return Rp2p1

    def get_ratio_AcAb_Aa(self):
        """ The function calculates the ratio of the difference between the b-point amplitude and c-point amplitude to the a-point amplitude.

            :return: (Ac-Ab)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ac = self.segment_d2[self.get_c()]
        return (Ac-Ab)/Aa

    def get_ratio_AdAb_Aa(self):
        """ The function calculates the ratio of the difference between the b-point amplitude and d-point amplitude to the a-point amplitude.

            :return: (Ad-Ab)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ad = self.segment_d2[self.get_d()]
        return (Ad-Ab)/Aa

    def getAGI(self):
        """ The function calculates the Aging Index.

            :return: (Ab-Ac-Ad-Ae)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ac = self.segment_d2[self.get_c()]
        Ad = self.segment_d2[self.get_d()]
        Ae = self.segment_d2[self.get_e()]
        return (Ab-Ac-Ad-Ae)/Aa

    def getAGImod(self):
        """ The function calculates the Modified Aging Index.

            :return: (Ab-Ac-Ad)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ac = self.segment_d2[self.get_c()]
        Ad = self.segment_d2[self.get_d()]
        return (Ab-Ac-Ad)/Aa

    def getAGIinf(self):
        """ The function calculates the Informal Aging Index.

            :return: (Ab-Ae)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ae = self.segment_d2[self.get_e()]
        return (Ab-Ae)/Aa

    def getAI(self):
        """ The function calculates the Augmentation Index, (PPG(Tp2) − PPG(Tp1))/Asp.

            :return: AI biomarker
        """
        AI = (self.segment[self.p2]-self.segment[self.p1])/self.peak_value
        return AI

    def getRIp1(self):
        """ The function calculates the Reflection Index of p1, Adp/(PPG(Tp1) − PPG(Tpi(0))).

            :return: RIp1 biomarker
        """
        RIp1 = self.getDiastolicPeak()/self.segment[self.p1]
        return RIp1

    def getRIp2(self):
        """ The function calculates the Reflection Index of p2, Adp/(PPG(Tp2) − PPG(Tpi(0))).

            :return: RIp2 biomarker
        """
        RIp2 = self.getDiastolicPeak()/self.segment[self.p2]
        return RIp2

    def getSC(self):
        """ The function calculates the Spring Constant, PPG"(Tsp)/((Asp-Au)/Asp).

            :return: SC biomarker
        """
        ddxSPA=self.segment_d2[(self.getTsp()*self.sample_rate).astype(int)]
        SPA=self.getSystolicPeak()
        MS=self.segment[self.u]
        SC = ddxSPA/((SPA-MS)/SPA)
        return SC

    def getIPAD(self):
        """ The function calculates the Inflection point area plus normalised d-point amplitude, AUCdia/AUCsys+Ad/Aa.

            :return: IPAD biomarker
        """
        IPAD = self.getAUCdia()/self.getAUCsys()+self.get_ratio_Ad_Aa()
        return IPAD
    def getRatioSW_DW(self, d):
        """ The function calculates the ratio of systolic and diastolic width at d percent of the pulse height

            :param d: the percentage chosen to calculate the width

            :return: ratio biomarker
        """
        sw_d = self.getSystolicWidth_d_percent(d)
        dw_d = self.getDiastolicWidth_d_percent(d)
        ratio = dw_d/sw_d
        return ratio

    def getPIR(self):
        """ The function calculates the ratio between the peak value and the right onset value

            :return: pir biomarker
        """
        pir = self.peak_value/self.onsets_values[1]
        return pir

    def getMS(self):
        """ The function calculates Maximum slope, PPG'(u)/(PPG(systolic peak) − PPG(systolic onset))

            :return: MS biomarker
        """
        MS = self.segment_d1[self.u]
        return MS

    def getUpslope(self):
        """ The function calculates Systolic Upslope between the left onset and the systolic peak.

            :return: Systolic Upslope
        """
        left_onset_time = self.onsets_times[0]*self.sample_rate
        left_onset_value = self.onsets_values[0]
        slope_numer = self.peak_value-left_onset_value
        slope_denom = self.peak_time*self.sample_rate - left_onset_time
        return slope_numer/slope_denom

    def getdiffVal(self):
        """ The function calculates the time between the left onset and the systolic peak.

            :return: left onset and systolic peak time
        """
        left_onset_value = self.onsets_values[0]
        diff = self.peak_value - left_onset_value
        return diff

    def getSTT(self):
        """ STT means slope transit time, which based on geometrical considerations of the PPG pulse wave to account for simultaneous.

            :return: STT biomarker
        """
        upslope = self.getUpslope()
        A = self.getdiffVal()
        return A/upslope

    def getSystolicPeakTime(self):
        """ Systolic Peak Time means the distance between the consecutive Systolic Peaks

        :return: Systolic Peak Times
         """
        return self.peak_time - self.onsets_times[0]

    def getSystolicPeakOutputCurve(self):
        """Peak time divided by systolic amplitude

        :return: sys_peak_time/sys_amplitude
        """
        sys_peak_time = self.getSystolicPeakTime()
        sys_amplitude = self.getSystolicPeak()
        return sys_peak_time/sys_amplitude

    def getAGIext(self):
        """ The function calculates the Extended Aging Index.

            :return: (Ab-Ac-Ad-Ae-Af)/Aa biomarker
        """
        Aa = self.segment_d2[self.get_a()]
        Ab = self.segment_d2[self.get_b()]
        Ac = self.segment_d2[self.get_c()]
        Ad = self.segment_d2[self.get_d()]
        Ae = self.segment_d2[self.get_e()]
        Af = self.segment_d2[self.get_e()]
        return (Ab-Ac-Ad-Ae-Af)/Aa

###########################################################################
############################ Get PPG biomarkers ###########################
###########################################################################
def get_biomarkers(s: pyPPG.PPG, fp: pyPPG.Fiducials, biomarkers_lst):
    """
    The function calculates the biomedical biomarkers of PPG signal.

    :param s: object of PPG signal
    :type s: pyPPG.PPG object
    :param fp: object of fiducial points
    :type fp: pyPPG.Fiducials object

    :return:
        - df: data frame with onsets, offset and peaks
        - df_biomarkers: data frame with PPG signal biomarkers
    """

    fs=s.fs
    ppg=s.filt_sig
    data = DotMap()

    df = pd.DataFrame()
    df_biomarkers = pd.DataFrame(columns=biomarkers_lst)
    peaks = fp.sp.values
    onsets = fp.on.values
    offsets = fp.off.values

    for i in range(len(onsets)):
        onset = onsets[i]
        offset = offsets[i]
        data.sig = ppg[int(onset):int(offset)]
        data.d1 = s.filt_d1[int(onset):int(offset)]
        data.d2 = s.filt_d2[int(onset):int(offset)]
        data.d3 = s.filt_d3[int(onset):int(offset)]
        peak = peaks[(peaks > onset) * (peaks < offset)]
        if len(peak) != 1:
            continue
        peak = peak[0]

        temp_fiducials = fp.get_row(i)

        peak_value = ppg[peak]
        peak_time = peak / fs
        onset_value = ppg[onset]
        onset_time = onset / fs

        if (peak_value - onset_value) == 0:
            continue

        offset_value = ppg[offset]
        offset_time = offset / fs

        idx_array = np.where(peaks == peak)
        idx = idx_array[0]
        onsets_values = np.array([onset_value, offset_value])
        onsets_times = np.array([onset_time, offset_time])
        if (idx + 1) < len(peaks):
            next_peak_value = ppg[peaks[idx + 1].astype('int64')][0]
            next_peak_time = peaks[idx + 1] / fs
            next_peak_time = next_peak_time[0]
            try:
                biomarkers_extractor = BmExctator(data, peak_value, peak_time, next_peak_value, next_peak_time, onsets_values, onsets_times, fs, biomarkers_lst,temp_fiducials)
                biomarkers_vec = biomarkers_extractor.get_biomarker_extract_func()
                lst = list(biomarkers_vec)
                df_biomarkers.loc[len(df_biomarkers.index)] = lst
                df = pd.concat({'onset': onset, 'offset': offset, 'peak': peak}, ignore_index=True)
            except:
                pass
        else:
            print("no more peaks")
    return df, df_biomarkers