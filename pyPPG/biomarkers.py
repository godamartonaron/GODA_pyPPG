import pyPPG

from pyPPG.ppg_bm.ppg_sig import get_ppg_sig
from pyPPG.ppg_bm.sig_ratios import get_sig_ratios
from pyPPG.ppg_bm.ppg_derivs import get_ppg_derivs
from pyPPG.ppg_bm.derivs_ratios import get_derivs_ratios
from pyPPG.ppg_bm.statistics import get_statistics

class BmCollection:

    ###########################################################################
    ######################## Initialization of Biomarkers #####################
    ###########################################################################
    def __init__(self, s: pyPPG.PPG, fp: pyPPG.Fiducials):
        """
        The purpose of the Biomarkers class is to calculate the PPG biomarkers.

        :param s: object of PPG signal
        :type s: pyPPG.PPG object
        :param fp: object of fiducial points
        :type fp: pyPPG.Fiducials object

        """

        self.s = s
        self.fp = fp

    ###########################################################################
    ############################ Get PPG Biomarkers ###########################
    ###########################################################################
    def get_biomarkers (self):
        """
        This function retrieves the list of biomarkers, computes their values, and calculates associated statistics.

        :return:
            - bm_defs: dictionary of biomarkers with name, definition and unit
            - bm_vals: dictionary of biomarkers with values
            - bm_stats: dictionary of biomarkers with statistics
        """

        s=self.s
        fp = self.fp

        ## Get Biomarkers
        bm_ppg_sig, def_ppg_sig = get_ppg_sig(s, fp)
        bm_sig_ratios, def_sig_ratios = get_sig_ratios(s, fp)
        bm_ppg_derivs, def_ppg_derivs = get_ppg_derivs(s, fp)
        bm_derivs_ratios, def_derivs_ratios = get_derivs_ratios(s, fp)

        bm_vals={'ppg_sig': bm_ppg_sig , 'sig_ratios': bm_sig_ratios, 'ppg_derivs': bm_ppg_derivs, 'derivs_ratios': bm_derivs_ratios}
        bm_defs = {'ppg_sig': def_ppg_sig, 'sig_ratios': def_sig_ratios, 'ppg_derivs': def_ppg_derivs, 'derivs_ratios': def_derivs_ratios}

        ## Get Statistics
        bm_stats = get_statistics(fp.sp, fp.on, bm_vals)

        ## Update index names
        BM_keys = bm_vals.keys()
        for key in BM_keys:
            bm_vals[key] = bm_vals[key].rename_axis('Index of pulse')
            bm_defs[key] = bm_defs[key].rename_axis('No. biomarkers')
            bm_stats[key] = bm_stats[key].rename_axis('Statistics')

        return bm_defs, bm_vals, bm_stats
