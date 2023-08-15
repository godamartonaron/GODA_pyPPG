import pyPPG

import pandas as pd
from pyPPG.ppg_bm.bm_extraction import get_biomarkers

###########################################################################
##################### Get Biomarkers of PPG Derivatives ###################
###########################################################################
def get_ppg_derivs(s: pyPPG.PPG, fp: pyPPG.Fiducials):
    """
    This function returns the biomarkers of PPG derivatives.

    :param s: object of PPG signal
    :type s: pyPPG.PPG object
    :param fp: object of fiducial points
    :type fp: pyPPG.Fiducials object

    :return:
        - biomarkers: dictionary of biomarkers of PPG derivatives
        - biomarkers_lst: list a biomarkers with name, definition and unit
    """

    biomarkers_lst = [
                    ["Tu",       "u-point time, the time between the pulse onset and u-point", "[s]"],
                    ["Tv",       "v-point time, the time between the pulse onset and v-point", "[s]"],
                    ["Tw",       "w-point time, the time between the pulse onset and w-point", "[s]"],
                    ["Ta",       "a-point time, the time between the pulse onset and a-point", "[s]"],
                    ["Tb",       "b-point time, the time between the pulse onset and b-point", "[s]"],
                    ["Tc",       "c-point time, the time between the pulse onset and c-point", "[s]"],
                    ["Td",       "d-point time, the time between the pulse onset and d-point", "[s]"],
                    ["Te",       "e-point time, the time between the pulse onset and e-point", "[s]"],
                    ["Tf",       "f-point time, the time between the pulse onset and f-point", "[s]"],
                    ["Tb-c",	 "b-c time, the time between the b-point and c-point", "[s]"],
                    ["Tb-d",	 "b-d time, the time between the b-point and d-point", "[s]"],
                    ["Tp1",	     "p1-point time, the time between the pulse onset and p1-point", "[s]"],
                    ["Tp2",      "p2-point time, the time between the pulse onset and p2-point", "[s]"],
                    ["Tp1-dp",   "p1-dia time, the time between the p1-point and diastolic peak", "[s]"],
                    ["Tp2-dp",   "p2-dia time, the time between the p2-point and diastolic peak", "[s]"],
    ]

    header = ['name', 'definition', 'unit']
    biomarkers_lst = pd.DataFrame(biomarkers_lst, columns=header)

    df, df_biomarkers = get_biomarkers(s, fp, biomarkers_lst.name)
    
    return df_biomarkers, biomarkers_lst
