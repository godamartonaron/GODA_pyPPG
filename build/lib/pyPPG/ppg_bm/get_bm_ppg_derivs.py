from pyPPG.ppg_bm.get_biomarkers import*

###########################################################################
##################### Get Biomarkers of PPG Derivatives ###################
###########################################################################
def get_bm_ppg_derivs(s, fiducials):
    """
    This function returns the biomarkers of PPG derivatives.

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
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.
        - PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
        - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
        - 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)

    :return: dictionary of biomarkers of PPG derivatives
    """

    biomarkers_lst = ["Tu",       # u-point time, the time between the pulse onset and u-point
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

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst)

    return df_biomarkers
