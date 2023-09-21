from pyPPG.pack_ppg._ErrorHandler import _check_shape_, WrongParameter
import pandas as pd

class PPG:
    '''
    This is class for the input PPG parameters.
    '''
    def __init__(self,s={}, check_ppg_len=True):
        """
        :param s: dictionary of the PPG signal:

            * s.start_sig: beginning of the signal in sample
            * s.end_sig: end of the signal in sample
            * s.v: a vector of PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.name: name of the record
            * s.v: 1-d array, a vector of raw PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.ppg: 1-d array, a vector of the filtered PPG values
            * s.vpg: 1-d array, a vector of the filtered PPG' values
            * s.apg: 1-d array, a vector of the filtered PPG" values
            * s.jpg: 1-d array, a vector of the filtered PPG'" values
            * s.filtering: a bool for filtering
            * s.correction: DataFrame where the key is the name of the fiducial points and the value is bool
        :type s: DotMap
        :param check_ppg_len: a bool for checking ppg length and sampling frequency
        :type check_ppg_len: bool

        """

        if s.fs <= 0:
            raise WrongParameter("Sampling frequency should be strictly positive")

        if check_ppg_len: _check_shape_(s.v, s.fs)

        s.check_ppg_len=check_ppg_len

        try:
            s.start_sig>0
        except:
            s.start_sig = 0

        try:
            s.end_sig>-1
        except:
            s.end_sig = -1

        # Initialise the correction for fiducial points
        if len(s.correction)<1:
            corr_on = ['on', 'dn', 'dp', 'v', 'w', 'f']
            corr_off = ['dn']
            correction=pd.DataFrame()
            correction.loc[0, corr_on] = True
            correction.loc[0, corr_off] = False
            s.correction=correction

        keys=s.keys()
        keys_list = list(keys)
        for i in keys_list:
            exec('self.'+i+' = s[i]')

    def get_s(self):
        """
        This function retrieves the dictionary of the PPG signal.

        :return: s: dictionary of the PPG signal
        """
        keys = self.__dict__.keys()
        keys_list = list(keys)
        s = {}
        for i in keys_list:
            s[i] = getattr(self, i, None)

        return pd.DataFrame(s)


class Fiducials:
        '''
        This is class for the PPG fiducial points.
        '''
        def __init__(self, fp: pd.DataFrame):
            """
            :param fiducials: DataFrame where the key is the name of the fiducial points and the value is the list of fiducial points PPG Fiducials Points

                * PPG signal (fp.on, fp.sp, fp.dn, fp.dp): List of pulse onset, systolic peak, dicrotic notch, diastolic peak
                * 1st derivative (fp.u, fp.v, fp.w): List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals
                * 2nd derivative (fp.a, fp.b, fp.c, fp.d, fp.e, fp.f): List of maximum and minimum points in 2nd derivitive between the onset to onset intervals
                * 3rd derivative (fp.p1, fp.p2): List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals
            :type fiducials: DataFrame

            """

            keys = fp.keys()
            keys_list = list(keys)
            for i in keys_list:
                exec('self.' + i + ' = fp[i]')

        def get_fp(self):
            """
            This function retrieves the struct of the fiducial points.

            :return: fp: DataFrame of the fiducial points
            """

            keys = self.__dict__.keys()
            keys_list = list(keys)
            fp={}
            for i in keys_list:
                fp[i]=getattr(self, i, None)

            return pd.DataFrame(fp)

        def get_row(self, row_index: int):
            """
            This function retrieves the specified row from the DataFrame of fiducial points.

            :param row_index: the index corresponding to the row in the fiducial points DataFrame
            :type row_index: int

            :return: the corresponding row in the fiducial points DataFrame
            """
            keys = self.__dict__.keys()
            keys_list = list(keys)
            fp_row={}
            for i in keys_list:
                fp_row[i]=getattr(self, i, None)[row_index]

            return pd.DataFrame(fp_row,[row_index])

class Biomarkers:
        '''
        This is class for the PPG biomarkers.
        '''
        def __init__(self, bm_defs={}, bm_vals={}, bm_stats={}):
            """
                This class constitutes a comprehensive dictionary encompassing biomarker definitions, values, and statistics. Each dictionary is organized into the subsequent subdirectories:
                    * ppg_sig: description for the PPG signal
                    * sig_ratios: description for the Signal ratios
                    * ppg_derivs: description for the PPG derivatives
                    * derivs_ratios: description for the Derivatives ratios

                :param bm_defs: dictionary with name, definition and unit of biomarkers in different categories:
                :type bm_defs: dict
                :param bm_vals: dictionary with values of biomarkers in different categories:
                :type bm_vals: dict
                :param bm_stats: data frame with summary of PPG features
                :type bm_stats: dict

            """

            for i in ['bm_defs', 'bm_vals', 'bm_stats']:
                exec('self.' + i + ' = ' +i)

        def get_bm(self):
            """
            This function retrieves the dictionary of the biomarkers.

            :return: bm: dictionary of the biomarkers
            """

            keys = self.__dict__.keys()
            keys_list = list(keys)
            bm={}
            for i in keys_list:
                bm[i]=getattr(self, i, None)

            return pd.DataFrame(bm)
