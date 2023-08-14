from pyPPG.pack_ppg._ErrorHandler import _check_shape_, WrongParameter
import pandas as pd

class PPG:
    '''
    This is class for the input PPG parameters.
    '''
    def __init__(self,s):
        """
        :param s: dictionary  of the PPG signal:

            * s.start: beginning of the signal in sample
            * s.end: end of the signal in sample
            * s.v: a vector of PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.name: name of the record
            * s.v: 1-d array, a vector of PPG values
            * s.fs: the sampling frequency of the PPG in Hz
            * s.filt_sig: 1-d array, a vector of the filtered PPG values
            * s.filt_d1: 1-d array, a vector of the filtered PPG' values
            * s.filt_d2: 1-d array, a vector of the filtered PPG" values
            * s.filt_d3: 1-d array, a vector of the filtered PPG'" values
        :type s: DotMap

        """

        if s.fs <= 0:
            raise WrongParameter("Sampling frequency should be strictly positive")
        _check_shape_(s.v, s.fs)

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
        def __init__(self, fp):
            """
            :param fiducials: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.

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

        def get_row(self, row_index):
            """
            This function retrieves the specified row from the DataFrame of fiducial points.

            :param row_index: the index corresponding to the row in the fiducial points DataFrame

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
        def __init__(self, bm_defs, bm_vals, bm_stats):
            """
                :param bm_defs: dictionary with name, definition and unit of biomarkers in different categories:

                    * def_ppg_sig: description of the PPG signal
                    * def_sig_ratios: description of the Signal ratios
                    * def_ppg_derivs: description of the PPG derivatives
                    * def_derivs_ratios: description of the Derivatives ratios
                :type bm_defs: dict
                :param bm_vals: dictionary with values of biomarkers in different categories:

                    * bm_ppg_sig: biomarkers of the PPG signal
                    * bm_sig_ratios: biomarkers of the Signal ratios
                    * bm_ppg_derivs: biomarkers of the PPG derivatives
                    * bm_derivs_ratios: biomarkers of the Derivatives ratios
                :type bm_vals: dict
                :param bm_stats: data frame with summary of PPG features

                    * def_ppg_sig: description of the PPG signal
                    * def_sig_ratios: description of the Signal ratios
                    * def_ppg_derivs: description of the PPG derivatives
                    * def_derivs_ratios: description of the Derivatives ratios
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
