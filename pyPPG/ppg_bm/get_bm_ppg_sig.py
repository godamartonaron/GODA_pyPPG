
from pyPPG.ppg_bm.biomarker_extractor import*

###########################################################################
####################### Get Biomarkers of PPG Signal ######################
###########################################################################
def get_bm_ppg_sig(s: pyPPG.PPG, fiducials: pd.DataFrame):
    """
    This function returns the biomarkers of PPG signal.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :param fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points PPG Fiducials Points.

        * PPG signal: List of pulse onset, systolic peak, dicrotic notch, diastolic peak
        * 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        * 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)
        * 3rd derivative: List of points of 1st maximum and minimum in 3rd derivitive between the onset to onset intervals (p1, p2)
    :type fiducials: DataFrame

    :return biomarkers: dictionary of biomarkers of PPG signal
    :return biomarkers_lst: list a biomarkers with name, definition and unit
    """

    biomarkers_lst = [
                    ["Tpi",   "Pulse interval, the time between the pulse onset and pulse offset", "[s]"],
                    ["Tpp",   "Peak-to-peak interval, the time between two consecutive systolic peaks", "[s]"],
                    ["Tsys",  "Systolic time, the time between the pulse onset and dicrotic notch", "[s]"],
                    ["Tdia",  "Diastolic time, the time is between the dicrotic notch and pulse offset", "[s]"],
                    ["Tsp",   "Systolic peak time, the time between the pulse onset and systolic peak", "[s]"],
                    ["Tdp",	  "Diastolic peak time, the time between the pulse onset and diastolic peak", "[s]"],
                    ["deltaT","Time delay, the time between the systolic peak and diastolic peak", "[s]"],
                    ["Tsw10", "Systolic width, the width at 10% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw25", "Systolic width, the width at 25% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw33", "Systolic width, the width at 33% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw50", "Systolic width, the width at 50% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw66", "Systolic width, the width at 66% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw75", "Systolic width, the width at 75% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tsw90", "Systolic width, the width at 90% of the systolic peak amplitude between the pulse onset and systolic peak", "[s]"],
                    ["Tdw10", "Diastolic width, the width at 10% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw25", "Diastolic width, the width at 25% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw33", "Diastolic width, the width at 33% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw50", "Diastolic width, the width at 50% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw66", "Diastolic width, the width at 66% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw75", "Diastolic width, the width at 75% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tdw90", "Diastolic width, the width at 90% of the systolic peak amplitude between the systolic peak and pulse offset", "[s]"],
                    ["Tpw10", "Pulse width, the sum of the systolic width and the diastolic width at 10%", "[s]"],
                    ["Tpw25", "Pulse width, the sum of the systolic width and the diastolic width at 25%", "[s]"],
                    ["Tpw33", "Pulse width, the sum of the systolic width and the diastolic width at 33%", "[s]"],
                    ["Tpw50", "Pulse width, the sum of the systolic width and the diastolic width at 50%", "[s]"],
                    ["Tpw66", "Pulse width, the sum of the systolic width and the diastolic width at 66%", "[s]"],
                    ["Tpw75", "Pulse width, the sum of the systolic width and the diastolic width at 75%", "[s]"],
                    ["Tpw90", "Pulse width, the sum of the systolic width and the diastolic width at 90%", "[s]"],
                    ["Asp",   "Systolic peak amplitude, the difference in amplitude between the pulse onset and systolic peak", "[nu]"],
                    ["Adn",   "Dicrotic notch amplitude, the difference in amplitude between the pulse onset and dicrotic notch", "[nu]"],
                    ["Adp",   "Diastolic peak amplitude, the difference in amplitude between the pulse onset and diastolic peak", "[nu]"],
                    ["Aoff",  "Pulse onset amplitude, the difference in amplitude between the pulse onset and pulse offset", "[nu]"],
                    ["AUCpi", "Area under pulse interval curve, the area under the pulse wave between pulse onset and pulse offset", "[nu]"],
                    ["AUCsys","Area under systolic curve, the area under the pulse wave between the pulse onset and the dicrotic notch", "[nu]"],
                    ["AUCdia","Area under diastolic curve, the area under the pulse wave between the dicrotic notch and pulse offset", "[nu]"],
    ]

    header = ['name', 'definition', 'unit']
    biomarkers_lst = pd.DataFrame(biomarkers_lst, columns=header)

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst.name)

    return df_biomarkers, biomarkers_lst
