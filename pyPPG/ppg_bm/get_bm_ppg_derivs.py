from pyPPG.ppg_bm.biomarker_extractor import*

###########################################################################
##################### Get Biomarkers of PPG Derivatives ###################
###########################################################################
def get_bm_ppg_derivs(s: pyPPG.PPG, fiducials: pd.DataFrame):
    """
    This function returns the biomarkers of PPG derivatives.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.

        * PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
        * 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        * 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
        * 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)
    :type fiducials: DataFrame

    :return biomarkers: dictionary of biomarkers of PPG derivatives
    :return biomarkers_lst: list a biomarkers with name, definition and unit
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

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst.name)
    
    return df_biomarkers, biomarkers_lst
