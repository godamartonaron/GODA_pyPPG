import copy
import pandas as pd

from pack_ppg._ErrorHandler import _check_shape_, WrongParameter
import numpy as np
from dotmap import DotMap
from scipy.signal import kaiserord, firwin, filtfilt, detrend, periodogram, lfilter, find_peaks, firls, resample
import matplotlib.pyplot as plt
import time
from scipy import signal


class FiducialPoints:
    def __init__(self, signal: np.array, fs: int, n_pools: int = 1):
        """
        The purpose of the FiducialPoints class is to calculate the fiducial points.

        :param signal: The PPG signal as a two-dimensional ndarray, when the first dimension is the len of the ecg, and the second is the number of leads.
        :param fs: The sampling frequency of the signal.
        :param peaks: The indexes of the R- points of the ECG signal – optional input.
        :param n_pools: The number of cores to use when calculating the fiducials.

        .. code-block:: python

            from pack_ppg.pyppg import FiducialPoints as Fp
            fp = Fp.FiducialPoints(f_ecg_rec, fs)

        """
        if fs <= 0:
            raise WrongParameter("Sampling frequency should be strictly positive")
        _check_shape_(signal, fs)

        self.signal = signal
        self.fs = fs
        self.peaks = []
        if n_pools is None:
            self.n_pools = 1
        else:
            self.n_pools = n_pools

###########################################################################
############################ Get Fiducials Points #########################
###########################################################################
def getFiducialsPoints(s):
    '''The function calculates the PPG Fiducials Points.
        - Original signal: List of pulse onset, pea and dicrotic notch
        - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
        - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values

    :return fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points.
    '''

    peak_detector='abp'

    drt0_fp=pd.DataFrame()
    peaks, onsets = abdp_beat_detector(s, peak_detector)
    dicroticnotch = getDicroticNotch(s, peaks, onsets)
    diastolicpeak = getDiastolicPeak(s, peaks, onsets, dicroticnotch)

    keys=('os', 'pk', 'dn','dp')
    dummy = np.empty(len(peaks))
    dummy.fill(np.NaN)
    n=0
    for temp_val in (onsets,peaks,dicroticnotch,diastolicpeak):
        drt0_fp[keys[n]] = dummy
        drt0_fp[keys[n]][0:len(temp_val)]=temp_val
        n=n+1

    drt1_fp = getFirstDerivitivePoints(s, onsets)
    drt2_fp = getSecondDerivitivePoints(s, onsets)
    drt3_fp = getThirdDerivitivePoints(s, onsets, drt2_fp)

    fiducials=pd.DataFrame()
    for temp_drt in (drt0_fp,drt1_fp,drt2_fp,drt3_fp):
        for key in list(temp_drt.keys()):
            fiducials[key] = dummy
            temp_val = temp_drt[key].values
            fiducials[key][0:len(temp_val)]=temp_val

    fiducials=fiducials.astype('Int64')
    return fiducials


###########################################################################
############################ PPG beat detector ############################
###########################################################################
def abdp_beat_detector(s, peak_detector):
    '''ABDP_BEAT_DETECTOR detects beats in a photoplethysmogram (PPG) signal
    using the improved 'Automatic Beat Detection' beat detector of Aboy M et al.

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values


    :return:
        - peaks: indices of detected pulse peaks
        - onsets: indices of detected pulse onsets

    Reference
    ---------
    Aboy M et al., An automatic beat detection algorithm for pressure signals.
    IEEE Trans Biomed Eng 2005; 52: 1662 - 70. <https://doi.org/10.1109/TBME.2005.855725>

    # Author
    Marton A. Goda: Faculty of Biomedical Engineering,
    Technion – Israel Institute of Technology, Haifa, Israel (October 2022)

    # Original Matlab implementation:
   Peter H. Charlton: King's College London (August 2017) – University of Cambridge (February 2022)
   <https://github.com/peterhcharlton/ppg-beats>

    # Changes from Charlton's implementation:
    I) Detect Maxima:
        1)  Systolic peak-to-peak distance is predicted by the heart rate estimate
            over the preceding 10 sec window.
        2)  The peak location is estimated by distances and prominences of the previous peaks.
    II) Find Onsets:
        1)  The onset is a local minimum, which is always calculated
            from the peak that follows it within a given time window
    III) Tidy of Peaks and Onsets:
        1)  There is a one-to-one correspondence between onsets and peaks
        2)  There are only onset and peak pairs
        3)  The distance between the onset and peak pairs can't be smaller than 30 ms
    '''

    # inputs
    x = copy.deepcopy(s.filt_sig)                          #signal
    fso=s.fs
    fs = 75
    x = resample(x, int(len(s.filt_sig)*(fs/fso)))
    up = setup_up_abdp_algorithm()                  #settings
    win_sec=10
    w = fs * win_sec                                #window length(number of samples)
    win_starts = np.array(list(range(0,len(x),round(0.8*w))))
    win_starts = win_starts[0:min(np.where([win_starts >= len(x) - w])[1])]
    win_starts = np.insert(win_starts,len(win_starts), len(x) + 1 - w)

    # before pre-processing
    hr_win=0  #the estimated systolic peak-to-peak distance, initially it is 0
    hr_win_v=[]
    px = DetectMaxima(x, 0, hr_win, peak_detector) # detect all maxima
    if len(px)==0:
        peaks = []
        return

    # detect peaks in windows
    all_p4 = []
    all_hr = np.empty(len(win_starts)-1)
    all_hr [:] = np.NaN
    hr_past = 0 # the actual heart rate
    hrvi = 0    # heart rate variability index

    for win_no in range(0,len(win_starts) - 1):
        curr_els = range(win_starts[win_no],win_starts[win_no] + w)
        curr_x = x[curr_els]

        y1 = Bandpass(curr_x, fs, 0.9 * up.fl_hz, 3 * up.fh_hz)     #Filter no.1
        hr = EstimateHeartRate(y1, fs, up, hr_past)  #Estimate HR from weakly filtered signal
        hr_past=hr
        all_hr[win_no] = hr

        if (peak_detector=='abp') and (hr>40):
            if win_no==0:
                p1 = DetectMaxima(y1, 0, hr_win, peak_detector)
                tr = np.percentile(np.diff(p1), 50)
                pks_diff = np.diff(p1)
                pks_diff = pks_diff[pks_diff>=tr]
                hrvi = np.std(pks_diff) / np.mean(pks_diff) * 5

            hr_win = fs / ((1 + hrvi) * 3)
            hr_win_v.append(hr_win)
        else:
            hr_win=0

        y2 = Bandpass(curr_x, fs, 0.9 * up.fl_hz, 2.5 * hr / 60)    # Filter no.2
        y2_deriv = EstimateDeriv(y2)    #Estimate derivative from highly filtered signal
        p2 = DetectMaxima(y2_deriv, up.deriv_threshold,hr_win, peak_detector) #Detect maxima in derivative
        y3 = Bandpass(curr_x, fs, 0.9 * up.fl_hz, 10 * hr / 60)
        p3 = DetectMaxima(y3, 50, hr_win, peak_detector)   #Detect maxima in moderately filtered signal
        p4 = find_pulse_peaks(p2, p3)
        p4 = np.unique(p4)

        if peak_detector=='abp':
            if len(p4)>round(win_sec/2):
                pks_diff = np.diff(p4)
                tr = np.percentile(pks_diff, 30)
                pks_diff = pks_diff[pks_diff >= tr]

                med_hr=np.median(all_hr[np.where(all_hr>0)])
                if ((med_hr*0.5<np.mean(pks_diff)) and (med_hr*1.5<np.mean(pks_diff))):
                    hrvi = np.std(pks_diff) / np.mean(pks_diff)*10

        all_p4 = np.concatenate((all_p4, win_starts[win_no] + p4), axis=None)

    all_p4=all_p4.astype(int)
    all_p4 = np.unique(all_p4)

    # IBT_0 = time.time()
    # if len(all_p4)>0:
    #     peaks, fn = IBICorrect(all_p4, px, np.median(all_hr), fs, up)
    #     peaks = np.unique(peaks)
    # else:
    #     peaks = all_p4
    # print('IBICorrect Time: ' + str(time.time() - IBT_0))

    peaks = (all_p4/fs*fso).astype(int)
    onsets, peaks = find_onsets(s.filt_sig, fso, up, peaks,60/np.median(all_hr)*fs)

    temp_i = np.where(np.diff(onsets) == 0)[0]
    if len(temp_i) > 0:
        peaks = np.delete(peaks, temp_i)
        onsets = np.delete(onsets, temp_i)

    temp_i = np.where((peaks - onsets) < fso / 30)[0]
    if len(temp_i) > 0:
        peaks = np.delete(peaks, temp_i)
        onsets = np.delete(onsets, temp_i)

    return peaks, onsets

###########################################################################
############################# Maximum detector ############################
###########################################################################
def DetectMaxima(sig, percentile,hr_win, peak_detector):
    #Table VI pseudocode
    """
    Detect Maxima function detects all peaks in the raw and also in the filtered signal to find.

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param percentile: in each signal partition, a rank filter detects the peaks above a given percentile
    :type percentile: int
    :param hr_win: window for adaptive the heart rate estimate
    :type hr_win: int

    :return: maximum peaks of signal, 1-d numpy array.

    """

    tr = np.percentile(sig, percentile)
    ld = len(sig)

    if peak_detector=='aby':

        s1,s2,s3 = sig[2:], sig[1:-1],sig[0:-2]
        m = 1 + np.array(np.where((s1 < s2) & (s3 < s2)))
        max_pks = m[sig[m] > tr]

    if peak_detector=='abp':
        s1,s2,s3 = sig[2:], sig[1:-1],sig[0:-2]

        max_loc = []
        min_loc = []
        max_pks=[]
        intensity_v = []
        if hr_win == 0:
            m = 1 + np.array(np.where((s1 < s2) & (s3 < s2)))
            max_pks = m[sig[m] > tr]
        else:
            max_loc = find_peaks(sig, distance=hr_win)[0]
            min_loc = find_peaks(-sig, distance=hr_win)[0]

            for i in range(0,len(max_loc)):
                values = abs(max_loc[i] - min_loc)
                min_v = min(values)
                min_i = np.where(min_v==values)[0][0]
                intensity_v.append(sig[max_loc[i]] - sig[min_loc[min_i]])

            # improvements:
            #   - adaptive threshold
            #   - probability density of maximum

            tr2 = np.mean(intensity_v)*0.25
            max_pks = find_peaks(sig+min(sig),prominence=tr2,distance=hr_win)[0]

    return max_pks

###########################################################################
############################ Bandpass filtering ###########################
###########################################################################
def Bandpass(sig, fs, lower_cutoff, upper_cutoff):
    """
    Bandpass filter function detects all peaks in the raw and also in the filtered signal to find.

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param fs: sampling frequency
    :type fs: int
    :param lower_cutoff: lower cutoff frequency
    :type lower_cutoff: float
    :param upper_cutoff: upper cutoff frequency
    :type upper_cutoff: float

    :return: bandpass filtered signal, 1-d numpy array.

    """

    # Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
    up = DotMap()
    up.paramSet.elim_vlf.Fpass = 1.3*lower_cutoff   #in Hz
    up.paramSet.elim_vlf.Fstop = 0.8*lower_cutoff   #in Hz
    up.paramSet.elim_vlf.Dpass = 0.05
    up.paramSet.elim_vlf.Dstop = 0.01

    # Filter characteristics: Eliminate VHFs (above frequency content of signals)
    up.paramSet.elim_vhf.Fpass = 1.2*upper_cutoff   #in Hz
    up.paramSet.elim_vhf.Fstop = 0.8*upper_cutoff   #in Hz
    up.paramSet.elim_vhf.Dpass = 0.05
    up.paramSet.elim_vhf.Dstop = 0.03

    # perform BPF
    s = DotMap()
    s.v = sig
    s.fs = fs

    b, a = signal.iirfilter(5, [2 * np.pi * lower_cutoff, 2 * np.pi * upper_cutoff], rs=60,
                            btype='band', analog=True, ftype='cheby2')

    bpf_sig = filtfilt(b, 1, s.v)

    return bpf_sig

###########################################################################
################### Filter the high frequency components  #################
###########################################################################
def elim_vlfs_abd(s, up, lower_cutoff):
    """
    This function filter the high frequency components.

    :param s: 1-d array, of shape (N,) where N is the length of the signal
    :param up: setup up parameters of the algorithm
    :type up: DotMap
    :param lower_cutoff: lower cutoff frequency
    :type lower_cutoff: float

    :return: high frequency filtered signal, 1-d numpy array.

    """

    ## Filter pre-processed signal to remove frequencies below resp
    # Adapted from RRest

    ## Eliminate nans
    s.v[np.isnan(s.v)] = np.mean(s.v[~np.isnan(s.v)])

    ##Make filter
    fc=lower_cutoff
    ripple=-20*np.log10(up.paramSet.elim_vlf.Dstop)
    width=abs(up.paramSet.elim_vlf.Fpass-up.paramSet.elim_vlf.Fstop)/(s.fs/2)
    [N,beta] = kaiserord(ripple,width)
    if N * 3 > len(s):
        N = round(N / 3)
    b = firwin(N, fc * 2 / s.fs, window=('kaiser', beta), scale=('True'))
    AMfilter = b#dfilt.dffir(b)
    #AMfilter = AMfilter[:50]
    ## Check frequency response
    # Gives a -3 dB cutoff at ? Hz, using:
    # freqz(AMfilter.Numerator)
    # norm_cutoff_freq = 0.0266;    % insert freq here from plot
    # cutoff_freq = norm_cutoff_freq*(fs/2);

    s_filt=DotMap()
    try:
        # s_filt.v = filtfilt(AMfilter.numerator, 1, s.v)
        s_filt.v = filtfilt(AMfilter, 1, s.v)
        s_filt.v = s.v-s_filt.v
    except:
        s_filt.v = s.v

    s_filt.fs = s.fs

    return s_filt

###########################################################################
################### Filter the low frequency components  ##################
###########################################################################
def elim_vhfs(s, up, upper_cutoff):
    """
    This function filter the high frequency components.

    :param s: 1-d array, of shape (N,) where N is the length of the signal
    :param up: setup up parameters of the algorithm
    :type up: DotMap
    :param upper_cutoff: upper cutoff frequency
    :type upper_cutoff: float

    :return: low frequency filtered signal, 1-d numpy array.

    """


    ## Filter signal to remove VHFs
    # Adapted from RRest
    s_filt = DotMap()

    ##Eliminate nans
    s.v[np.isnan(s.v)] = np.mean(s.v[~np.isnan(s.v)])

    ##Check to see if sampling freq is at least twice the freq of interest
    if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1:
        #then the fs is too low to perform this filtering
        s_filt.v = s.v
        return

    ## Create filter
    # parameters for the low-pass filter to be used
    # flag  = 'scale';

    fc = upper_cutoff
    ripple = -20 * np.log10(up.paramSet.elim_vhf.Dstop)
    width = abs(up.paramSet.elim_vhf.Fpass - up.paramSet.elim_vhf.Fstop) / (s.fs / 2)
    [N, beta] = kaiserord(ripple, width)
    if N * 3 > len(s):
        N = round(N / 3)
    #b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag)
    b = firwin(N, fc * 2 / s.fs, window=('kaiser', beta), scale=('True'))
    AMfilter = b#dfilt.dffir(b)
    #AMfilter = AMfilter[:50]

    ## Check frequency response
    # Gives a -3 dB cutoff at cutoff_freq Hz, using:
    # freqz(AMfilter.Numerator)
    # norm_cutoff_freq = 0.3355;    % insert freq here from plot
    # cutoff_freq = norm_cutoff_freq*(s.fs/2);

    ## Remove VHFs
    s_dt=detrend(s.v)
    s_filt.v = filtfilt(AMfilter, 1, s_dt)

    return s_filt

###########################################################################
########################### Heart Rate estimation #########################
###########################################################################
def EstimateHeartRate(sig, fs, up, hr_past):
    """
    Heart Rate Estimation function estimate the heart rate according to the previous heart rate in given time window

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param fs: sampling frequency
    :type fs: int
    :param up: setup up parameters of the algorithm
    :type up: DotMap
    :param hr_past: the average heart rate in the past in given time window
    :type hr_past: int

    :return: estimated heart rate, 1-d numpy array.

    """

    # Estimate PSD
    blackman_window = np.blackman(len(sig))
    f, pxx = periodogram(sig,fs, blackman_window)
    ph = pxx
    fh = f

    # Extract HR
    if (hr_past / 60 < up.fl_hz) | (hr_past / 60 > up.fh_hz):
        rel_els = np.where((fh >= up.fl_hz) & (fh <= up.fh_hz))
    else:
        rel_els = np.where((fh >= hr_past / 60 * 0.5) & (fh <= hr_past / 60 * 1.4))

    rel_p = ph[rel_els]
    rel_f = fh[rel_els]
    max_el = np.where(rel_p==max(rel_p))
    hr = rel_f[max_el]*60

    return hr

###########################################################################
############# Estimate derivative from highly filtered signal #############
###########################################################################
def EstimateDeriv(sig):
    """
    Derivative Estimation function estimate derivative from highly filtered signal based on the
    General least-squares smoothing and differentiation by the convolution (Savitzky Golay) method

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :return: derivative, 1-d numpy array.

    """

    #Savitzky Golay
    deriv_no = 1
    win_size = 5
    deriv = savitzky_golay_abd(sig, deriv_no, win_size)

    return deriv


def savitzky_golay_abd(sig, deriv_no, win_size):
    """
    This function estimate the Savitzky Golay derivative from highly filtered signal

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param deriv_no: number of derivative
    :type deriv_no: int
    :param win_size: size of window
    :type win_size: int

    :return: Savitzky Golay derivative, 1-d numpy array.

    """

    ##assign coefficients
    # From: https: // en.wikipedia.org / wiki / Savitzky % E2 % 80 % 93 Golay_filter  # Tables_of_selected_convolution_coefficients
    # which are calculated from: A., Gorry(1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method".Analytical Chemistry. 62(6): 570?3. doi: 10.1021 / ac00205a007.

    if deriv_no==0:
        #smoothing
        if win_size == 5:
            coeffs = [-3, 12, 17, 12, -3]
            norm_factor = 35
        elif win_size == 7:
            coeffs = [-2, 3, 6, 7, 6, 3, -2]
            norm_factor = 21
        elif win_size == 9:
            coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21]
            norm_factor = 231
        else:
            print('Can''t do this window size')
    elif deriv_no==1:
        # first derivative
        if win_size == 5:
            coeffs = range(-2,3)
            norm_factor = 10
        elif win_size == 7:
            coeffs = range(-3,4)
            norm_factor = 28
        elif win_size == 9:
            coeffs = range(-4,5)
            norm_factor = 60
        else:
            print('Can''t do this window size')
    elif deriv_no == 2:
        # second derivative
        if win_size == 5:
            coeffs = [2, -1, -2, -1, 2]
            norm_factor = 7
        elif win_size == 7:
            coeffs = [5, 0, -3, -4, -3, 0, 5]
            norm_factor = 42
        elif win_size == 9:
            coeffs = [28, 7, -8, -17, -20, -17, -8, 7, 28]
            norm_factor = 462
        else:
            print('Can''t do this window size')
    elif deriv_no == 3:
        # third derivative
        if win_size == 5:
            coeffs = [-1, 2, 0, -2, 1]
            norm_factor = 2
        elif win_size == 7:
            coeffs = [-1, 1, 1, 0, -1, -1, 1]
            norm_factor = 6
        elif win_size == 9:
            coeffs = [-14, 7, 13, 9, 0, -9, -13, -7, 14]
            norm_factor = 198
        else:
            print('Can''t do this window size')
    elif deriv_no == 4:
        # fourth derivative
        if win_size == 7:
            coeffs = [3, -7, 1, 6, 1, -7, 3]
            norm_factor = 11
        elif win_size == 9:
            coeffs = [14, -21, -11, 9, 18, 9, -11, -21, 14]
            norm_factor = 143
        else:
            print('Can''t do this window size')
    else:
        print('Can''t do this order of derivative')


    if deriv_no % 2 == 1:
        coeffs = -np.array(coeffs)

    A = [1, 0]
    filtered_sig = lfilter(coeffs, A, sig)
    # filtered_sig = filtfilt(coeffs, A, sig)
    s = len(sig)
    half_win_size = np.floor(win_size * 0.5)
    zero_pad=filtered_sig[win_size] * np.ones(int(half_win_size))
    sig_in=filtered_sig[win_size-1:s]
    sig_end=filtered_sig[s-1] * np.ones(int(half_win_size))
    deriv = [*zero_pad,*sig_in,*sig_end]
    deriv = deriv / np.array(norm_factor)

    return deriv

###########################################################################
############################# Pulse detection #############################
###########################################################################
def find_pulse_peaks(p2,p3):
    """
    Pulse detection function detect the pulse peaks according to the peaks of 1st and 2nd derivatives
    General least-squares smoothing and differentiation by the convolution (Savitzky Golay) method

    :param p2: 1-d array, peaks of the 1st derivatives
    :param p3: 1-d array, peaks of the 2nd derivatives
    :return: pulse peaks, 1-d numpy array.

    """

    p4 = np.empty(len(p2))
    p4[:] = np.NaN
    for k in range(0,len(p2)):
        rel_el = np.where(p3>p2[k])
        if np.any(rel_el) and ~np.isnan(rel_el[0][0]):
            p4[k] = p3[rel_el[0][0]]

    p4 = p4[np.where(~np.isnan(p4))]
    p4 = p4.astype(int)
    return p4

###########################################################################
####################### Correct peaks' location error #####################
###########################################################################
def  IBICorrect(p, m, hr, fs, up):
    """
    This function corrects the peaks' location error

    :param p: systolic peaks of the PPG signal
    :type p: 1-d numpy array
    :param m: all maxima of the PPG signal
    :type m: 1-d numpy array
    :param hr: heart rate
    :type hr: 1-d numpy array
    :param fs: sampling frequency
    :type fs: int
    :param up: setup up parameters of the algorithm
    :type up: DotMap

    :return: onsets, 1-d numpy array.

    """

    #Correct peaks' location error due to pre-processing
    pc = np.empty(len(p))
    pc[:] = np.NaN
    pc1=[]
    for k in range(0,len(p)):
        temp_pk=abs(m - p[k])
        rel_el = np.where(temp_pk==min(temp_pk))
        pc1=[*pc1,*m[rel_el]]

    # Correct false positives
    # identify FPs
    d = np.diff(pc1)/fs    # interbeat intervals in secs
    fp = find_reduced_IBIs(d, hr, up)
    # remove FPs
    pc = np.array(pc1)[fp]

    # Correct false negatives
    d = np.diff(pc)/fs    # interbeat intervals in secs
    fn = find_prolonged_IBIs(d, hr, up)

    pc = pc1

    return pc, fn

def find_reduced_IBIs(IBIs, med_hr, up):
    IBI_thresh = up.lower_hr_thresh_prop*60/med_hr
    fp = IBIs < IBI_thresh
    fp = [*np.where(fp == 0)[0].astype(int)]
    return fp

def find_prolonged_IBIs(IBIs, med_hr, up):
    IBI_thresh = up.upper_hr_thresh_prop*60/med_hr
    fn = IBIs > IBI_thresh
    # fn = [*np.where(fn == 0)[0].astype(int)]
    return fn

###########################################################################
####################### Setup up the beat detector ########################
###########################################################################
def setup_up_abdp_algorithm():
    """
    This function setups the filter parameters of the algorithm

    :return: filter parameters of the algorithm, DotMap.

    """
    # plausible HR limits
    up=DotMap()
    up.fl = 30               #lower bound for HR
    up.fh = 200              #upper bound for HR
    up.fl_hz = up.fl/60
    up.fh_hz = up.fh/60

    # Thresholds
    up.deriv_threshold = 75          #originally 90
    up.upper_hr_thresh_prop = 2.25   #originally 1.75
    up.lower_hr_thresh_prop = 0.5    #originally 0.75

    # Other parameters
    up.win_size = 10    #in secs

    return up

###########################################################################
############################## Find PPG onsets ############################
###########################################################################
def find_onsets(sig,fs,up,peaks,med_hr):
    """
    This function finds the onsets of PPG sigal

    :param sig: 1-d array, of shape (N,) where N is the length of the signal
    :param fs: sampling frequency
    :type fs: int
    :param up: setup up parameters of the algorithm
    :type up: DotMap
    :param peaks: peaks of the signal
    :type peaks: 1-d array

    :return: onsets, 1-d numpy array.

    """

    Y1=Bandpass(sig, fs, 0.9*up.fl_hz, 3*up.fh_hz)
    temp_oi0=find_peaks(-Y1,distance=med_hr*0.3)[0]

    null_indexes = np.where(temp_oi0<peaks[0])
    if len(null_indexes[0])!=0:
        if len(null_indexes[0])==1:
            onsets = [null_indexes[0][0]]
        else:
            onsets = [null_indexes[0][-1]]
    else:
        onsets = [peaks[0]-round(fs/50)]

    i=1
    while i < len(peaks):
        min_SUT=fs*0.12     # minimum Systolic Upslope Time 120 ms
        min_DT=fs*0.3       # minimum Diastolic Time 300 ms

        before_peak=temp_oi0 <peaks[i]
        after_last_onset=temp_oi0 > onsets[i - 1]
        SUT_time=peaks[i]-temp_oi0>min_SUT
        DT_time = temp_oi0-peaks[i-1]  > min_DT
        temp_oi1 = temp_oi0[np.where(before_peak * after_last_onset*SUT_time*DT_time)]
        if len(temp_oi1)>0:
            if len(temp_oi1) == 1:
                onsets.append(temp_oi1[0])
            else:
                onsets.append(temp_oi1[-1])
            i=i+1
        else:
            peaks = np.delete(peaks, i)

    return onsets,peaks

###########################################################################
########################## Detect dicrotic notch ##########################
###########################################################################
def getDicroticNotch (s, peaks, onsets):
    """
    Dicrotic Notch function estimate the location of dicrotic notch in between the systolic and diastolic peak

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param peaks: 1-d array, peaks of the signal
    :param onsets: 1-d array, onsets of the signal

    :return: location of dicrotic notches, 1-d numpy array.
    """

    ## The 2nd derivative and Hamming low pass filter is calculated.
    dxx = np.diff(np.diff(s.filt_sig))
    fs = s.fs

    # Make filter
    Fn = fs / 2                                 # Nyquist Frequency
    FcU = 20                                    # Cut off Frequency: 20 Hz
    FcD = FcU + 5                               # Transition Frequency: 5 Hz

    n = 21                                      # Filter order
    f = [0, (FcU / Fn), (FcD / Fn), 1]          # Frequency band edges
    a = [1, 1, 0, 0]                            # Amplitudes
    b = firls(n, f, a)

    lp_filt_sig = filtfilt(b, 1,  dxx)    # Low pass filtered signal with 20 cut off Frequency and 5 Hz Transition width

    ## The weighting is calculated and applied to each beat individually
    def t_wmax(i, peaks,onsets):
        if i < 3:
            HR = np.mean(np.diff(peaks))/fs
            t_wmax = -0.1 * HR + 0.45
        else:
            t_wmax = np.mean(peaks[i - 3:i]-onsets[i - 3:i])/fs
        return t_wmax

    dic_not=[]
    for i in range(0,len(onsets)-1):
        nth_beat = lp_filt_sig[onsets[i]:onsets[i + 1]]

        i_Pmax=peaks[i]-onsets[i]
        t_Pmax=(peaks[i]-onsets[i])/fs
        t=np.linspace(0,len(nth_beat)-1,len(nth_beat))/fs
        T_beat=(len(nth_beat)-1)/fs
        tau=(t-t_Pmax)/(T_beat-t_Pmax)
        tau[0:i_Pmax] = 0
        beta=5

        if len(peaks)>1:
            t_w=t_wmax(i, peaks, onsets)
        else:
            t_w=np.NaN

        if t_w!=T_beat:
            tau_wmax=(t_w-t_Pmax)/(T_beat-t_Pmax)
        else:
            tau_wmax=0.9

        alfa=(beta*tau_wmax-2*tau_wmax+1)/(1-tau_wmax)
        if (alfa > 4.5) or (alfa < 1.5):
            HR = np.mean(np.diff(peaks))/fs
            t_w = -0.1 * HR + 0.45
            tau_wmax = (t_w - t_Pmax) / (T_beat - t_Pmax)
            alfa = (beta * tau_wmax - 2 * tau_wmax + 1) / (1 - tau_wmax)

        ## Calculate the Dicrotic Notch for each heart cycle using the weighted window
        if alfa>1:
            w = tau ** (alfa - 1) * (1 - tau) ** (beta - 1)
        else:
            w = tau * (1 - tau) ** (beta - 1)

        pp=w*nth_beat
        pp = pp[np.where(~np.isnan(pp))]
        max_pp_v = np.max(pp)
        max_pp_i=np.where(pp==max_pp_v)[0][0]
        ## NOTE!! Shifting with 5 sample. FIX IT!
        dic_not.append(max_pp_i+onsets[i]+5)

    return dic_not

###########################################################################
########################## Detect diastolic peak ##########################
###########################################################################
def getDiastolicPeak(s, peaks, onsets, dicroticnotch):
    """
    Dicrotic Notch function estimate the location of dicrotic notch in between the systolic and diastolic peak

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param peaks: 1-d array, peaks of the signal
    :param onsets: 1-d array, onsets of the signal
    :param dicroticnotches: 1-d array, onsets of the signal

    :return diastolicpeak: location of dicrotic notches, 1-d numpy array.
    """

    diastolicpeak=dicroticnotch

    return diastolicpeak

###########################################################################
####################### Get First Derivitive Points #######################
###########################################################################
def getFirstDerivitivePoints(s, onsets):
    """Calculate first derivitive points u and v from the PPG' signal

    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param onsets: 1-d array, onsets of the signal

    :return
        - u: The greatest maximum peak between the left systolic onset and the systolic peak on PPG'
        - v: The lowest minimum pits between the systolic peak and the right systolic onset on PPG'
        - w: The first maximum peak after v on PPG'
    """

    dx = s.filt_d1
    u, v, w = [], [], []
    for i in range(0,len(onsets)-1):
        segment = dx[onsets[i]:onsets[i + 1]]

        # u fiducial point
        max_loc = np.argmax(dx[onsets[i]:onsets[i + 1]])+onsets[i]
        u.append(max_loc)

        # v fiducial point
        min_loc = np.argmin(dx[u[-1]:onsets[i+1]])+u[-1]-1
        v.append(min_loc)

        # w fiducial point
        temp_segment=dx[v[-1]:onsets[i+1]]
        max_locs, _ = find_peaks(temp_segment)
        if max_locs.size>0:
            max_w = max_locs[0]+v[-1]-1
        else:
            max_w = np.NaN
        w.append(max_w)

    drt1_fp = pd.DataFrame({"u":[], "v":[], "w":[]})
    drt1_fp.u, drt1_fp.v, drt1_fp.w = u, v, w
    return drt1_fp

###########################################################################
####################### Get Second Derivitive Points ######################
###########################################################################
def getSecondDerivitivePoints(s, onsets):
    """Calculate Second derivitive points a, b, c, d, e, and f from the PPG" signal
    :param s: a struct of PPG signal:
        - s.v: a vector of PPG values
        - s.fs: the sampling frequency of the PPG in Hz
        - s.filt_sig: a vector of PPG values
        - s.filt_d1: a vector of PPG values
        - s.filt_d2: a vector of PPG values
        - s.filt_d3: a vector of PPG values
    :param onsets: 1-d array, onsets of the signal

    :return
        - a: The first maximum peak between left systolic onset and systolic peak on PPG"
        - b: The first minimum pits after a on PPG"
        - c: The greatest maximum peak between b and e, or if no maximum peak is present then the inflection point on PPG"
        - d: The lowest minimum pits between c and e, or if no minimum pits is present then the inflection point on PPG"
        - e: The greatest maximum peak between the systolic peak and  the right systolic onset
        - f: The first minimum pits after e and before 0.8 CP (cardiac period) on PPG"
    """

    sig = s.filt_sig
    ddx = s.filt_d2
    dddx = s.filt_d3

    a, b, c, d, e, f = [], [], [], [], [], []
    for i in range(0,len(onsets)-1):

        # a fiducial point
        temp_pk=np.argmax(sig[onsets[i]:onsets[i + 1]])+onsets[i]-1
        temp_segment=ddx[onsets[i]:temp_pk]
        max_locs, _ = find_peaks(temp_segment)
        max_loc = max_locs[np.argmax(temp_segment[max_locs])]
        max_a = max_loc + onsets[i] - 1
        a.append(max_a)

        # b fiducial point
        temp_segment=ddx[a[-1]:onsets[i+1]]
        min_locs, _ = find_peaks(-temp_segment)
        min_b = min_locs[0]+a[-1]-1
        b.append(min_b)

        # e fiducial point
        e_lower_bound = b[-1]
        upper_bound_coeff = 0.6
        e_upper_bound = ((onsets[i + 1] - onsets[i]) * upper_bound_coeff + onsets[i]).astype(int)
        temp_segment=ddx[e_lower_bound:e_upper_bound]
        max_locs, _ = find_peaks(temp_segment)
        if max_locs.size==0:
            temp_segment=ddx[b[-1]:onsets[i + 1]]
            max_locs, _ = find_peaks(temp_segment)

        max_loc = max_locs[np.argmax(temp_segment[max_locs])]
        max_e = max_loc + e_lower_bound - 1

        e.append(max_e)

        # c fiducial point
        temp_segment=ddx[b[-1]:e[-1]]
        max_locs, _ = find_peaks(temp_segment)
        if max_locs.size>0:
            max_loc = max_locs[0]
        else:
            temp_segment=dddx[b[-1]:e[-1]]
            min_locs, _ = find_peaks(-temp_segment)

            if min_locs.size > 0:
                max_loc = min_locs[np.argmin(temp_segment[min_locs])]
            else:
                max_locs, _ = find_peaks(temp_segment)
                max_loc = max_locs[0]

        max_c = max_loc + b[-1] - 1
        c.append(max_c)

        # d fiducial point
        temp_segment = ddx[c[-1]:e[-1]]
        min_locs, _ = find_peaks(-temp_segment)
        if min_locs.size > 0:
            min_loc = min_locs[np.argmin(temp_segment[min_locs])]
            min_d = min_loc + c[-1] - 1
        else:
            min_d = max_c

        d.append(min_d)

        # f fiducial point
        temp_segment = ddx[e[-1]:onsets[i + 1]]
        min_locs, _ = find_peaks(-temp_segment)
        if (min_locs.size > 0) and (min_locs[0]<len(sig)*0.8):
            min_loc = min_locs[0]
        else:
            min_loc = []

        min_f = min_loc + e[-1] - 1
        f.append(min_f)

    drt2_fp = pd.DataFrame({"a":[], "b":[], "c":[],"d":[], "e":[], "f":[]})
    drt2_fp.a, drt2_fp.b, drt2_fp.c, drt2_fp.d, drt2_fp.e, drt2_fp.f = a, b, c, d, e, f
    return drt2_fp

def getThirdDerivitivePoints(s, onsets,drt2_fp):
    """Calculate third derivitive points p1 and p2 from the PPG'" signal
        :param s: a struct of PPG signal:
            - s.v: a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.filt_sig: a vector of PPG values
            - s.filt_d1: a vector of PPG values
            - s.filt_d2: a vector of PPG values
            - s.filt_d3: a vector of PPG values
        :param onsets: 1-d array, onsets of the signal

        :return
        - p1: The first local maximum of PPG'"  after b.
        - p2: Identify a candidate p2 at the last local minimum of PPG'"  before d (unless c=d, in which case take the
        first local minimum of PPG'"  after d). If there is a local maximum of PPG'" between this candidate and dic then
        use this instead.

    """
    dddx = s.filt_d3

    p1, p2 = [], []
    for i in range(0,len(onsets)-1):
        # p1 fiducial point
        ref_b=drt2_fp.b[np.squeeze(np.where(np.logical_and(drt2_fp.b>onsets[i],drt2_fp.b<onsets[i+1])))]
        if ref_b.size==0:
            ref_b=onsets[i]

        temp_segment=dddx[ref_b:onsets[i+1]]
        max_locs, _ = find_peaks(temp_segment)
        max_loc = max_locs[0]
        max_p1 = max_loc + ref_b - 1
        p1.append(max_p1)

        # p2 fiducial point
        ref_start = p1[-1]
        if ref_start.size>0:
            ref_c = drt2_fp.c[np.squeeze(np.where(np.logical_and(drt2_fp.c>onsets[i],drt2_fp.c<onsets[i+1])))]
            ref_d = drt2_fp.d[np.squeeze(np.where(np.logical_and(drt2_fp.d>onsets[i],drt2_fp.d<onsets[i+1])))]

            if ref_d>ref_c:
                ref_end=ref_d
                min_ind=-1
            elif ref_c.size>0:
                ref_start=ref_c
                ref_end = onsets[i+1]
                min_ind = 0
            elif drt2_fp.e.size>0:
                ref_end = onsets[i+1]
                min_ind = 0

            temp_segment=dddx[ref_start:ref_end]
            min_locs, _ = find_peaks(-temp_segment)
            # min_p2 = min_locs[min_ind] + ref_start - 1
            if min_locs.size>0:
                min_p2 = min_locs[min_ind]+ref_start-1
            else:
                min_p2 = ref_end

            p2.append(min_p2)
        else:
            p2.append(ref_start)

    drt3_fp = pd.DataFrame({"p1": [], "p2": []})
    drt3_fp.p1, drt3_fp.p2 = p1, p2

    return drt3_fp
