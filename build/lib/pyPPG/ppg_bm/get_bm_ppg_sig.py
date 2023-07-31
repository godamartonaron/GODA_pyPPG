from pyPPG.ppg_bm.get_biomarkers import*

###########################################################################
####################### Get Biomarkers of PPG Signal ######################
###########################################################################
def get_bm_ppg_sig(s, fiducials):
    """
    This function returns the biomarkers of PPG signal.

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

    :return: dictionary of biomarkers of PPG signal
    """

    biomarkers_lst = ["Tpi",   # Pulse Interval, the time between the pulse onset and pulse offset
                    "Tpp",   # Peak-to-Peak Interval, the time between two consecutive systolic peaks
                    "Tsys",	 # Systolic Time, the time between the pulse onset and dicrotic notch
                    "Tdia",  # Diastolic Time, the time is between the dicrotic notch and pulse offset
                    "Tsp",   # Systolic Peak Time, the time between the pulse onset and systolic peak
                    "Tdp",	 # Diastolic Peak Time, the time between the pulse onset and diastolic peak
                    "deltaT",# Time Delay, the time between the systolic peak and diastolic peak
                    "Tsw10", # Systolic Width, the width at 10% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw25", # Systolic Width, the width at 25% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw33", # Systolic Width, the width at 33% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw50", # Systolic Width, the width at 50% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw66", # Systolic Width, the width at 66% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw75", # Systolic Width, the width at 75% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tsw90", # Systolic Width, the width at 90% of the Systolic Peak Amplitude between the pulse onset and systolic peak
                    "Tdw10", # Diastolic Width, the width at 10% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw25", # Diastolic Width, the width at 25% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw33", # Diastolic Width, the width at 33% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw50", # Diastolic Width, the width at 50% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw66", # Diastolic Width, the width at 66% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw75", # Diastolic Width, the width at 75% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tdw90", # Diastolic Width, the width at 90% of the Systolic Peak Amplitude between the systolic peak and pulse offset
                    "Tpw10", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 10%
                    "Tpw25", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 25%
                    "Tpw33", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 33%
                    "Tpw50", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 50%
                    "Tpw66", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 66%
                    "Tpw75", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 75%
                    "Tpw90", # Pulse Width, the sum of the Systolic Width and the Diastolic Width at 90%
                    "Asp",   # Systolic Peak Amplitude, the difference in amplitude between the pulse onset and systolic peak
                    "Adn",   # Dicrotic Notch Amplitude, the difference in amplitude between the pulse onset and dicrotic notch
                    "Adp",   # Diastolic Peak Amplitude, the difference in amplitude between the pulse onset and diastolic peak
                    "Aoff",  # Pulse Onset Amplitude, the difference in amplitude between the pulse onset and pulse offset
                    "AUCpi", # Area Under Pulse Interval Curve, the area under the pulse wave between pulse onset and pulse offset
                    "AUCsys",# Area Under Systolic Curve, the area under the pulse wave between the pulse onset and the dicrotic notch
                    "AUCdia",# Area Under Diastolic Curve, the area under the pulse wave between the dicrotic notch and pulse offset
                    ]
    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst)

    return df_biomarkers
