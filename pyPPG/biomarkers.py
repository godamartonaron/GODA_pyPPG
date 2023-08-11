from ppg_bm.get_bm_ppg_sig import*
from ppg_bm.get_bm_sig_ratios import*
from ppg_bm.get_bm_ppg_derivs import*
from ppg_bm.get_bm_derivs_ratios import*

class Biomarkers:

    ###########################################################################
    ######################## Initialization of Biomarkers #####################
    ###########################################################################
    def __init__(self, s: pyPPG.PPG, fiducials: pd.DataFrame):
        """
        The purpose of the Biomarkers class is to calculate the ppg biomarkers.

        :param s: a struct of PPG signal
        :type s: pyPPG.PPG object
        :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.

            * PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
            * 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
            * 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
            * 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)
        :type fiducials: DataFrame

        """

        self.s = s
        self.fiducials = fiducials

    ###########################################################################
    ############################ Get PPG Biomarkers ###########################
    ###########################################################################
    def getBiomarkers (self):
        """
        This function returns the PPG biomarkers.

        :return: bm_vals: dictionary with values of biomarkers in different categories:

                    * bm_ppg_sig: biomarkers of the PPG signal
                    * bm_sig_ratios: biomarkers of the Signal ratios
                    * bm_ppg_derivs: biomarkers of the PPG derivatives
                    * bm_derivs_ratios: biomarkers of the Derivatives ratios

                bm_defs: dictionary with name, definition and unit of biomarkers in different categories:

                    * def_ppg_sig: description of the PPG signal
                    * def_sig_ratios: description of the Signal ratios
                    * def_ppg_derivs: description of the PPG derivatives
                    * def_derivs_ratios: description of the Derivatives ratios
        """

        s=self.s
        fiducials = self.fiducials

        bm_ppg_sig, def_ppg_sig = get_bm_ppg_sig(s,fiducials)
        bm_sig_ratios, def_sig_ratios = get_bm_sig_ratios(s, fiducials)
        bm_ppg_derivs, def_ppg_derivs = get_bm_ppg_derivs(s,fiducials)
        bm_derivs_ratios, def_derivs_ratios = get_bm_derivs_ratios(s,fiducials)

        bm_vals={'ppg_sig': bm_ppg_sig , 'sig_ratios': bm_sig_ratios, 'ppg_derivs': bm_ppg_derivs, 'derivs_ratios': bm_derivs_ratios}
        bm_defs = {'ppg_sig': def_ppg_sig, 'sig_ratios': def_sig_ratios, 'ppg_derivs': def_ppg_derivs, 'derivs_ratios': def_derivs_ratios}

        return bm_vals, bm_defs
