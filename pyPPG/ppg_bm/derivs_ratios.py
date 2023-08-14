from pyPPG.ppg_bm.bm_extractor import*

###########################################################################
####################### Get Ratios of PPG Derivatives #####################
###########################################################################
def get_derivs_ratios(s: pyPPG.PPG, fp: pyPPG.Fiducials):
    """
    This function returns the biomarkers of Derivatives ratios.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :param fp: a struct of fiducial points
    :type fp: pyPPG.Fiducials object

    :return:
        - df_biomarkers: dictionary of biomarkers of Derivatives ratios
        - biomarkers_lst: list a biomarkers with name, definition and unit
    """

    biomarkers_lst = [
                    ["Tu/Tpi",       "Ratio of the u-point time vs. the pulse interval", "[%]"],
                    ["Tv/Tpi",       "Ratio of the v-point time vs. the pulse interval", "[%]"],
                    ["Tw/Tpi",       "Ratio of the w-point time vs. the pulse interval", "[%]"],
                    ["Ta/Tpi",       "Ratio of the a-point time vs. the pulse interval", "[%]"],
                    ["Tb/Tpi",       "Ratio of the b-point time vs. the pulse interval", "[%]"],
                    ["Tc/Tpi",       "Ratio of the c-point time vs. the pulse interval", "[%]"],
                    ["Td/Tpi",       "Ratio of the d-point time vs. the pulse interval", "[%]"],
                    ["Te/Tpi",       "Ratio of the e-point time vs. the pulse interval", "[%]"],
                    ["Tf/Tpi",       "Ratio of the f-point time vs. the pulse interval", "[%]"],
                    ["(Tu-Ta)/Tpi",  "Ratio of the difference between the u-point time and a-point time vs. the pulse interval", "[%]"],
                    ["(Tv-Tb)/Tpi",  "Ratio of the difference between the v-point time and b-point time vs. the pulse interval", "[%]"],
                    ["Au/Asp",       "Ratio of the u-point amplitude vs. the systolic peak amplitude", "[%]"],
                    ["Av/Au",        "Ratio of the v-point amplitude vs. the u-point amplitude", "[%]"],
                    ["Aw/Au",        "Ratio of the w-point amplitude vs. the u-point amplitude", "[%]"],
                    ["Ab/Aa",        "Ratio of the b-point amplitude vs. the a-point amplitude", "[%]"],
                    ["Ac/Aa",        "Ratio of the c-point amplitude vs. the a-point amplitude", "[%]"],
                    ["Ad/Aa",        "Ratio of the d-point amplitude vs. the a-point amplitude", "[%]"],
                    ["Ae/Aa",        "Ratio of the e-point amplitude vs. the a-point amplitude", "[%]"],
                    ["Af/Aa",        "Ratio of the f-point amplitude vs. the a-point amplitude", "[%]"],
                    ["Ap2/Ap1",      "Ratio of the p2-point amplitude vs. the p1-point amplitude", "[%]"],
                    ["(Ac-Ab)/Aa",   "Ratio of the difference between the b-point amplitude and c-point amplitude vs. the a-point amplitude", "[%]"],
                    ["(Ad-Ab)/Aa",   "Ratio of the difference between the b-point amplitude and d-point amplitude vs. the a-point amplitude", "[%]"],
                    ["AGI",          "Aging Index, (Ab-Ac-Ad-Ae)/Aa", "[%]"],
                    ["AGImod",       "Modified aging index, (Ab-Ac-Ad)/Aa", "[%]"],
                    ["AGIinf",       "Informal aging index, (Ab-Ae)/Aa", "[%]"],
                    ["AI",           "Augmentation index, (PPG(Tp2)-PPG(Tp1))/Asp", "[%]"],
                    ["RIp1",         "Reflection index of p1, Adp/(PPG(Tp1)-PPG(Tpi(0)))", "[%]"],
                    ["RIp2",         "Reflection index of p2, Adp/(PPG(p2)-PPG(Tpi(0)))", "[%]"],
                    ["SC",           "Spring constant, PPG''(Tsp)/((Asp-Au)/Asp)", "[nu]"],
                    ["IPAD",         "Inflection point area plus normalised d-point amplitude, AUCdia/AUCsys+Ad/Aa", "[nu]"],
    ]

    header = ['name', 'definition', 'unit']
    biomarkers_lst = pd.DataFrame(biomarkers_lst, columns=header)

    df, df_biomarkers = get_biomarkers(s, fp, biomarkers_lst.name)

    return df_biomarkers, biomarkers_lst