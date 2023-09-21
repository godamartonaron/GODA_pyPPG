import numpy as np
import pandas as pd
import scipy.stats
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

###########################################################################
###################### Statistics of PPG Biomarkers #######################
###########################################################################

def get_statistics(peaks: pd.Series, onsets: pd.Series, ppg_biomarkers: dict):
    """
    The function compares the different biomedical features of PPG signal.

    :param peaks: 1-d array, peaks of the signal
    :type peaks: Series
    :param onsets: 1-d array, onsets of the signal
    :type onsets: Series
    :param ppg_biomarkers: dictionary of the PPG biomarkers
    :type ppg_biomarkers: dict

    :return: df_windows: data frame with summary of PPG features
    """

    BM_keys = list(ppg_biomarkers.keys())
    ppg_statistics={}
    for j in range (0, len(BM_keys)):
        df_features=ppg_biomarkers[BM_keys[j]]

        df_stat = pd.DataFrame()
        for k in df_features.keys():
            data = df_features[k].values
            df_tempstat = {}

            try: df_tempstat['mean'] = np.mean(data)
            except: pass

            try: df_tempstat['median'] = np.median(data)
            except: pass

            try: df_tempstat['std'] = np.std(data)
            except: pass

            try: percentile_25 = np.percentile(data, 25)
            except: pass

            try: percentile_75 = np.percentile(data, 75)
            except: pass

            try: df_tempstat['percentile_25'] = percentile_25
            except: pass

            try: df_tempstat['percentile_75'] = percentile_75
            except: pass

            try: df_tempstat['iqr'] = percentile_75 - percentile_25
            except: pass

            try: df_tempstat['skew'] = scipy.stats.skew(data)
            except: pass

            try: df_tempstat['kurtosis'] = scipy.stats.kurtosis(data)
            except: pass

            try: df_tempstat['mad'] = np.mean(np.absolute(data - np.mean(data)))
            except: pass

            try:df_stat[k]=df_tempstat.values()
            except: pass

        df_stat.index = ['mean', 'median', 'std', 'percentile_25', 'percentile_75', 'iqr', 'skew', 'kurtosis', 'mad']
        ppg_statistics[BM_keys[j]]=df_stat

    return ppg_statistics