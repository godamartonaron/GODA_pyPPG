import numpy as np
import pandas as pd
import scipy.stats

###########################################################################
########################## Summary of PPG signal ##########################
###########################################################################

def Summary(ppg, peaks, onsets, fs):
    """
    The function compares the different biomedical features of PPG signal.

    :param ppg: 1-d array, of shape (N,) where N is the length of the signal
    :param peaks: 1-d array, peaks of the signal
    :param onsets: 1-d array, onsets of the signal
    :param fs: sampling frequency
    :type fs: int

    :return df_windows: data frame with summary of PPG features
    """

    df = pd.DataFrame({"onset": onsets[:-1], "offset": onsets[1:], "peak": peaks[:-1]})
    #df, df_features = get_features(ppg, peaks, onsets, fs, features_lst)
    df_data = pd.concat([df], axis=1)
    #   remove outliers
    df = df_data
    df['duration'] = df['offset'] - df['onset']
    df.head(20)
    df = df.sort_values(by='duration')

    ppg_sum = {}
    ppg_sum['mean'] = np.mean(ppg)
    ppg_sum['median'] = np.median(ppg)
    ppg_sum['std'] = np.std(ppg)
    percentile_25 = np.percentile(ppg, 25)
    percentile_75 = np.percentile(ppg, 75)
    ppg_sum['percentile_25'] = percentile_25
    ppg_sum['percentile_75'] = percentile_75
    ppg_sum['iqr'] = percentile_75 - percentile_25
    ppg_sum['skew'] = scipy.stats.skew(ppg)
    ppg_sum['kurtosis'] = scipy.stats.kurtosis(ppg)
    ppg_sum['mad'] = np.mean(np.absolute(ppg - np.mean(ppg)))
    ppg_sum['n_peaks'] = len(peaks)
    ppg_sum = pd.DataFrame(ppg_sum.values(),ppg_sum.keys())

    return ppg_sum