from pyPPG.ppg_bm.biomarker_extractor import*

###########################################################################
######################### Get Ratios of PPG Signal ########################
###########################################################################
def get_bm_sig_ratios(s: pyPPG.PPG, fiducials: pd.DataFrame):
    """
    This function returns the biomarkers of Signal ratios.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object

    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.

        * PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
        * 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        * 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
        * 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)

    :return biomarkers: dictionary of biomarkers of Signal ratios
    :return biomarkers_lst: list a biomarkers with name, definition and unit
    """

    biomarkers_lst = [
                    ["IPR",          "Instantaneous pulse rate, 60 / Tpi", "[%]"],
                    ["Tsys/Tdia",    "Ratio of the systolic time vs. the diastolic time", "[%]"],
                    ["Tpw25/Tpi",    "Ratio of the pulse width at 25% of the systolic peak amplitude vs. the pulse interval", "[%]"],
                    ["Tpw50/Tpi",    "Ratio of the pulse width at 50% of the systolic peak amplitude vs. the pulse interval", "[%]"],
                    ["Tpw75/Tpi",    "Ratio of the pulse width at 75% of the systolic peak amplitude vs. the pulse interval", "[%]"],
                    ["Tpw25/Tsp",    "Ratio of the pulse width at 25% of the systolic peak amplitude vs. the systolic peak time", "[%]"],
                    ["Tpw50/Tsp",    "Ratio of the pulse width at 50% of the systolic peak amplitude vs. the systolic peak time", "[%]"],
                    ["Tpw75/Tsp",    "Ratio of the pulse width at 75% of the systolic peak amplitude vs. the systolic peak time", "[%]"],
                    ["Tdw10/Tsw10",  "Ratio of the diastolic width vs. the systolic width at 10% width", "[%]"],
                    ["Tdw25/Tsw25",  "Ratio of the diastolic width vs. the systolic width at 25% width", "[%]"],
                    ["Tdw33/Tsw33",  "Ratio of the diastolic width vs. the systolic width at 33% width", "[%]"],
                    ["Tdw50/Tsw50",  "Ratio of the diastolic width vs. the systolic width at 50% width", "[%]"],
                    ["Tdw66/Tsw66",  "Ratio of the diastolic width vs. the systolic width at 66% width", "[%]"],
                    ["Tdw75/Tsw75",  "Ratio of the diastolic width vs. the systolic width at 75% width", "[%]"],
                    ["Tdw90/Tsw90",  "Ratio of the diastolic width vs. the systolic width at 90% width", "[%]"],
                    ["Tsp/Tpi",      "Ratio of the systolic peak time vs. the pulse interval", "[%]"],
                    ["Asp/Aoff",     "Ratio of the systolic peak amplitude vs. the pulse offset amplitude", "[%]"],
                    ["Adp/Asp",      "Reflection index, the ratio of the diastolic peak amplitude vs. the systolic peak amplitude", "[%]"],
                    ["IPA",          "Inflection point area, the ratio of the area under diastolic curve vs. the area under systolic curve", "[nu]"],
                    ["Tsp/Asp",      "Ratio of the systolic peak time vs. the systolic peak amplitude", "[nu]"],
                    ["Asp/deltaT",   "Stiffness index, the ratio of the systolic peak amplitude vs. the time delay", "[nu]"],
                    ["Asp/(Tpi-Tsp)","Ratio of the systolic peak amplitude vs. the difference between the pulse interval and systolic peak time ", "[nu]"],
    ]

    header = ['name', 'definition', 'unit']
    biomarkers_lst = pd.DataFrame(biomarkers_lst, columns=header)

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst.name)

    return df_biomarkers, biomarkers_lst
