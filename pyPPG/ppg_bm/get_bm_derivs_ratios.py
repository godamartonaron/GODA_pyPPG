from pyPPG.ppg_bm.get_biomarkers import*

###########################################################################
####################### Get Ratios of PPG Derivatives #####################
###########################################################################
def get_bm_derivs_ratios(s, fiducials):
    """
    This function returns the biomarkers of Derivatives ratios.

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

    :return biomarkers: dictionary of biomarkers of Derivatives ratios
    """

    biomarkers_lst = ["Tu/Tpi",       # The ratio of the u-point time to the Pulse Interval
                    "Tv/Tpi",       # The ratio of the v-point time to the Pulse Interval
                    "Tw/Tpi",       # The ratio of the w-point time to the Pulse Interval
                    "Ta/Tpi",       # The ratio of the a-point time to the Pulse Interval
                    "Tb/Tpi",       # The ratio of the b-point time to the Pulse Interval
                    "Tc/Tpi",       # The ratio of the c-point time to the Pulse Interval
                    "Td/Tpi",       # The ratio of the d-point time to the Pulse Interval
                    "Te/Tpi",       # The ratio of the e-point time to the Pulse Interval
                    "Tf/Tpi",       # The ratio of the f-point time to the Pulse Interval
                    "(Tu-Ta)/Tpi",  # The ratio of the difference between the u-point time and a-point time to the Pulse Interval
                    "(Tv-Tb)/Tpi",  # The ratio of the difference between the v-point time and b-point time to the Pulse Interval
                    "Au/Asp",       # The ratio of the u-point amplitude to the Systolic Peak Amplitude
                    "Av/Au",        # The ratio of the v-point amplitude to the u-point amplitude
                    "Aw/Au",        # The ratio of the w-point amplitude to the u-point amplitude
                    "Ab/Aa",        # The ratio of the b-point amplitude to the a-point amplitude
                    "Ac/Aa",        # The ratio of the c-point amplitude to the a-point amplitude
                    "Ad/Aa",        # The ratio of the d-point amplitude to the a-point amplitude
                    "Ae/Aa",        # The ratio of the e-point amplitude to the a-point amplitude
                    "Af/Aa",        # The ratio of the f-point amplitude to the a-point amplitude
                    "Ap2/Ap1",      # The ratio of the p2-point amplitude to the p1-point amplitude
                    "(Ac-Ab)/Aa",   # The ratio of the difference between the b-point amplitude and c-point amplitude to the a-point amplitude
                    "(Ad-Ab)/Aa",   # The ratio of the difference between the b-point amplitude and d-point amplitude to the a-point amplitude
                    "AGI",          # Aging Index, (Ab-Ac-Ad-Ae)/Aa
                    "AGImod",       # Modified Aging Index, (Ab-Ac-Ad)/Aa
                    "AGIinf",       # Informal Aging Index, (Ab-Ae)/Aa
                    "AI",           # Augmentation Index, (PPG(Tp2) − PPG(Tp1))/Asp
                    "RIp1",         # Reflection Index of p1, Adp/(PPG(Tp1) − PPG(Tpi(0)))
                    "RIp2",         # Reflection Index of p2, Adp/(PPG(p2) − PPG(Tpi(0)))
                    "SC",           # Spring Constant, PPG"(Tsp)/((Asp-Au)/Asp)
                    "IPAD",         # Inflection point area plus normalised d-point amplitude, AUCdia/AUCsys+Ad/Aa
                     ]

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst)

    return df_biomarkers