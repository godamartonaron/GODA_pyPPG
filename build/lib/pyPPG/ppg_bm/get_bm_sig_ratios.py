from pyPPG.ppg_bm.get_biomarkers import*

###########################################################################
######################### Get Ratios of PPG Signal ########################
###########################################################################
def get_bm_sig_ratios(s, fiducials):
    """
    This function returns the biomarkers of Signal ratios.

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

    :return biomarkers: dictionary of biomarkers of Signal ratios
    """

    biomarkers_lst = ["IPR",          # Instantaneous Pulse Rate, 60 / Tpi
                    "Tsys/Tdia",    # The ratio of the Systolic Time to the Diastolic Time
                    "Tpw25/Tpi",    # The ratio of the Pulse Width at 25% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw50/Tpi",    # The ratio of the Pulse Width at 50% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw75/Tpi",    # The ratio of the Pulse Width at 75% of the Systolic Peak Amplitude to the Pulse Interval
                    "Tpw25/Tsp",    # The ratio of the Pulse Width at 25% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tpw50/Tsp",    # The ratio of the Pulse Width at 50% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tpw75/Tsp",    # The ratio of the Pulse Width at 75% of the Systolic Peak Amplitude to the Systolic Peak Time
                    "Tdw10/Tsw10",  # The ratio of the Diastolic Width to the Systolic Width at 10% width
                    "Tdw25/Tsw25",  # The ratio of the Diastolic Width to the Systolic Width at 25% width
                    "Tdw33/Tsw33",  # The ratio of the Diastolic Width to the Systolic Width at 33% width
                    "Tdw50/Tsw50",  # The ratio of the Diastolic Width to the Systolic Width at 50% width
                    "Tdw66/Tsw66",  # The ratio of the Diastolic Width to the Systolic Width at 66% width
                    "Tdw75/Tsw75",  # The ratio of the Diastolic Width to the Systolic Width at 75% width
                    "Tdw90/Tsw90",  # The ratio of the Diastolic Width to the Systolic Width at 90% width
                    "Tsp/Tpi",      # The ratio of the Systolic Peak Time to the Pulse Interval
                    "Asp/Aoff",     # The ratio of the Systolic Peak Amplitude to the Pulse Offset Amplitude
                    "Adp/Asp",      # Reflection Index, the ratio of the Diastolic Peak Amplitude to the Systolic Peak Amplitude
                    "IPA",          # Inflection Point Area, the ratio of the Area Under Diastolic Curve to the Area Under Systolic Curve
                    "Tsp/Asp",      # The ratio of the Systolic Peak Time to the Systolic Peak Amplitude
                    "Asp/deltaT",   # Stiffness Index, the ratio of the Systolic Peak Amplitude to the Time Delay
                    "Asp/(Tpi-Tsp)",# The ratio of the Systolic Peak Amplitude to the difference between the Pulse Interval and Systolic Peak Time
    ]

    df, df_biomarkers = get_biomarkers(s, fiducials, biomarkers_lst)

    return df_biomarkers
