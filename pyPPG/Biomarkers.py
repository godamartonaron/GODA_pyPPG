from pack_ppg._ErrorHandler import _check_shape_, WrongParameter
from ppg_bm.get_bm_ppg_sig import*
from ppg_bm.get_bm_sig_ratios import*
from ppg_bm.get_bm_ppg_derivs import*
from ppg_bm.get_bm_derivs_ratios import*

class Biomarkers:

    ###########################################################################
    ######################## Initialization of Biomarkers #####################
    ###########################################################################
    def __init__(self, s: DotMap, fiducials: pd.DataFrame):
        """
        The purpose of the Biomarkers class is to calculate the ppg biomarkers.

        :param s: a struct of PPG signal:
            - s.v: a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.name: name of the record
            - s.v: 1-d array, a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.filt_sig: 1-d array, a vector of the filtered PPG values
            - s.filt_d1: 1-d array, a vector of the filtered PPG' values
            - s.filt_d2: 1-d array, a vector of the filtered PPG" values
            - s.filt_d3: 1-d array, a vector of the filtered PPG'" values
        :type s: DotMap
        :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.
            - PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
            - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
            - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
            - 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)
        :type fiducials: DataFrame

        """
        if s.fs <= 0:
            raise WrongParameter("Sampling frequency should be strictly positive")
        _check_shape_(s.v, s.fs)

        self.s = s
        self.fiducials = fiducials

    ###########################################################################
    ############################ Get PPG Biomarkers ###########################
    ###########################################################################
    def getBiomarkers (self):
        """
        This function returns the PPG biomarkers.

        :return: bm_vals: dictionary of biomarkers in different categories:
            - PPG signal
            - Signal ratios
            - PPG derivatives
            - Derivatives ratios
        :return: bm_defs: dictionary of biomarkers with name, definition and unit:
            - PPG signal
            - Signal ratios
            - PPG derivatives
            - Derivatives ratios
        """

        s=self.s
        fiducials = self.fiducials

        bm_ppg_sig, lst_ppg_sig = get_bm_ppg_sig(s,fiducials)
        bm_sig_ratios, lst_sig_ratios = get_bm_sig_ratios(s, fiducials)
        bm_ppg_derivs, lst_ppg_derivs = get_bm_ppg_derivs(s,fiducials)
        bm_derivs_ratios, lst_derivs_ratios = get_bm_derivs_ratios(s,fiducials)

        bm_vals={'ppg_sig': bm_ppg_sig , 'sig_ratios': bm_sig_ratios, 'ppg_derivs': bm_ppg_derivs, 'derivs_ratios': bm_derivs_ratios}
        bm_defs = {'ppg_sig': lst_ppg_sig, 'sig_ratios': lst_sig_ratios, 'ppg_derivs': lst_ppg_derivs, 'derivs_ratios': lst_derivs_ratios}

        return bm_vals, bm_defs
