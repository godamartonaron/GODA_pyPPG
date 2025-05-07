import pyPPG

import copy
import pandas as pd
import numpy as np
from dotmap import DotMap
from scipy.signal import kaiserord, firwin, filtfilt, detrend, periodogram, lfilter, find_peaks, firls, resample
from scipy import signal

class FpCollection:

    ###########################################################################
    ###################### Initialization of Fiducial Points ##################
    ###########################################################################
    def __init__(self, s: pyPPG.PPG):
        """
        The purpose of the FiducialPoints class is to calculate the fiducial points.

        :param s: object of PPG signal
        :type s: pyPPG.PPG object

        """

        keys=s.__dict__.keys()
        keys_list = list(keys)
        for i in keys_list:
            exec('self.'+i+' = s.'+i)

    ###########################################################################
    ############################ Get Fiducial Points ##########################
    ###########################################################################
    def get_fiducials(self, s: pyPPG.PPG):
        '''This function calculates the PPG Fiducial Points.
            - Original signal: List of systolic peak, pulse onset, dicrotic notch, and diastolic peak
            - 1st derivative: List of points of 1st maximum and minimum in 1st derivitive between the onset to onset intervals (u,v)
            - 2nd derivative: List of maximum and minimum points in 2nd derivitive between the onset to onset intervals (a, b, c, d, e)

        :param s: object of PPG signal
        :type s: pyPPG.PPG object

        :return: fiducial points: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points
        '''

        # 'ABD' refers the original Aboy peak detectorm and 'PPGdet' refers the improved version.
        peak_detector='PPGdet'

        # Extract Fiducial Points
        ppg_fp=pd.DataFrame()
        peaks, onsets = self.get_peak_onset(peak_detector)
        dicroticnotch = self.get_dicrotic_notch(peaks, onsets)

        vpg_fp = self.get_vpg_fiducials(onsets)
        apg_fp = self.get_apg_fiducials(onsets, peaks)
        jpg_fp = self.get_jpg_fiducials(onsets, apg_fp)

        diastolicpeak = self.get_diastolic_peak(onsets, dicroticnotch, apg_fp.e)

        # Merge Fiducial Points
        keys=('on', 'sp', 'dn','dp')
        dummy = np.empty(len(peaks))
        dummy.fill(np.NaN)
        n=0
        for temp_val in (onsets,peaks,dicroticnotch,diastolicpeak):
            ppg_fp[keys[n]] = dummy
#            ppg_fp[keys[n]][0:len(temp_val)]=temp_val
            ppg_fp.loc[0:len(temp_val)-1, keys[n]] = temp_val
            n=n+1

        fiducials=pd.DataFrame()
        for temp_sig in (ppg_fp,vpg_fp,apg_fp,jpg_fp):
            for key in list(temp_sig.keys()):
                fiducials[key] = dummy
                temp_val = temp_sig[key].values
#                fiducials[key][0:len(temp_val)]=temp_val
                fiducials.loc[:len(temp_val) - 1, key] = temp_val

        # Correct Fiducial Points
        fiducials=self.correct_fiducials(fiducials, s.correction)

        fiducials=fiducials.astype('Int64')

        # Extract pulse offsets
        offsets = copy.deepcopy(fiducials.on[1:])
        offsets.index = offsets.index - 1

        # Update index name
        fiducials = fiducials.drop(len(fiducials)-1)
        fiducials = fiducials.rename_axis('Index of pulse')

        # Add pulse offsets
        fiducials.insert(4, 'off', offsets)

        return fiducials


    ###########################################################################
    ############################ PPG beat detector ############################
    ###########################################################################
    def get_peak_onset(self, peak_detector='PPGdet'):
        '''PPGdet detects beats in a photoplethysmogram (PPG) signal
        using the improved 'Automatic Beat Detection' of Aboy M et al.

        :param peak_detector: type of peak detector
        :type peak_detector: str

        :return:
            - peaks, 1-d array: indices of detected systolic peaks
            - onsets, 1-d array: indices of detected pulse onsets

        Reference
        ---------
        Aboy M et al., An automatic beat detection algorithm for pressure signals.
        IEEE Trans Biomed Eng 2005; 52: 1662 - 70. <https://doi.org/10.1109/TBME.2005.855725>

        Author:
        Marton A. Goda – Faculty of Biomedical Engineering,
        Technion – Israel Institute of Technology, Haifa, Israel (August 2022)

        Original Matlab implementation:
        Peter H. Charlton – King's College London (August 2017) – University of Cambridge (February 2022)
        <https://github.com/peterhcharlton/ppg-beats>

        Changes from Charlton's implementation:
            1) Detect Maxima:
                *  Systolic peak-to-peak distance is predicted by the heart rate estimate over the preceding 10 sec window.
                *  The peak location is estimated by distances and prominences of the previous peaks.
            2) Find Onsets:
                *  The onset is a local minimum, which is always calculated from the peak that follows it within a given time window
            3) Tidy of Peaks and Onsets:
                *  There is a one-to-one correspondence between onsets and peaks
                *  There are only onset and peak pairs
                *  The distance between the onset and peak pairs can't be smaller than 30 ms
            4) Correct Peaks and Onsets:
                * The Peaks must be the highest amplitude between two consecutive pulse onsets, if not, then these are corrected
                * After the correction of Peaks, the Onsets are recalculated
        '''

        # inputs
        x = copy.deepcopy(self.ppg)                    #signal
        fso=self.fs
        fs = 75
        x = resample(x, int(len(self.ppg)*(fs/fso)))
        up = self.set_beat_detection()                 #settings
        win_sec=10
        w = fs * win_sec                                    #window length(number of samples)
        win_starts = np.array(list(range(0,len(x),round(0.8*w))))
        win_starts = win_starts[0:min(np.where([win_starts >= len(x) - w])[1])]
        win_starts = np.insert(win_starts,len(win_starts), len(x) + 1 - w)

        # before pre-processing
        hr_win=0  #the estimated systolic peak-to-peak distance, initially it is 0
        hr_win_v=[]
        px = self.detect_maxima(x, 0, hr_win, peak_detector) # detect all maxima
        if len(px)==0:
            peaks = []
            onsets = []
            return peaks, onsets

        # detect peaks in windows
        all_p4 = []
        all_hr = np.empty(len(win_starts)-1)
        all_hr [:] = np.NaN
        hr_past = 0 # the actual heart rate
        hrvi = 0    # heart rate variability index

        for win_no in range(0,len(win_starts) - 1):
            curr_els = range(win_starts[win_no],win_starts[win_no] + w)
            curr_x = x[curr_els]

            y1 = self.def_bandpass(curr_x, fs, 0.9 * up.fl_hz, 3 * up.fh_hz)   # Filter no.1
            hr = self.estimate_HR(y1, fs, up, hr_past)               # Estimate HR from weakly filtered signal
            hr_past=hr
            all_hr[win_no] = hr

            if (peak_detector=='PPGdet') and (hr>40):
                if win_no==0:
                    p1 = self.detect_maxima(y1, 0, hr_win, peak_detector)
                    tr = np.percentile(np.diff(p1), 50)
                    pks_diff = np.diff(p1)
                    pks_diff = pks_diff[pks_diff>=tr]
                    hrvi = np.std(pks_diff) / np.mean(pks_diff) * 5

                hr_win = fs / ((1 + hrvi) * 3)
                hr_win_v.append(hr_win)
            else:
                hr_win=0

            y2 = self.def_bandpass(curr_x, fs, 0.9 * up.fl_hz, 2.5 * hr / 60)           # Filter no. 2
            y2_deriv = self.estimate_deriv(y2)                                          # Estimate derivative from highly filtered signal
            p2 = self.detect_maxima(y2_deriv, up.deriv_threshold,hr_win, peak_detector) # Detect maxima in derivative
            y3 = self.def_bandpass(curr_x, fs, 0.9 * up.fl_hz, 10 * hr / 60)
            p3 = self.detect_maxima(y3, 50, hr_win, peak_detector)                      # Detect maxima in moderately filtered signal
            p4 = self.find_pulse_peaks(p2, p3)
            p4 = np.unique(p4)

            if peak_detector=='PPGdet':
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

        peaks, fn = self.correct_IBI(all_p4, px, np.median(all_hr), fs, up)

        peaks = (all_p4/fs*fso).astype(int)
        onsets, peaks = self.find_onsets(self.ppg, fso, up, peaks,60/np.median(all_hr)*fs)

        # Correct Peaks
        for i in range(0, len(peaks) - 1):
            max_loc = np.argmax(self.ppg[onsets[i]:onsets[i + 1]]) + onsets[i]
            if peaks[i] != max_loc:
                peaks[i] = max_loc

        # Correct onsets
        onsets, peaks = self.find_onsets(self.ppg, fso, up, peaks, 60 / np.median(all_hr) * fs)

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
    def detect_maxima(self, sig: np.array, percentile: int ,hr_win: int, peak_detector: str):
        #Table VI pseudocode
        """
        Detect Maxima function detects all peaks in the raw and also in the filtered signal to find.

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array
        :param percentile: in each signal partition, a rank filter detects the peaks above a given percentile
        :type percentile: int
        :param hr_win: window for adaptive the heart rate estimate
        :type hr_win: int
        :param peak_detector: type of peak detector
        :type peak_detector: str

        :return: maximum peaks of signal, 1-d array.

        """

        tr = np.percentile(sig, percentile)

        if peak_detector=='ABD':

            s1,s2,s3 = sig[2:], sig[1:-1],sig[0:-2]
            m = 1 + np.array(np.where((s1 < s2) & (s3 < s2)))
            max_pks = m[sig[m] > tr]

        if peak_detector=='PPGdet':
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

                # possible improvements:
                #   - using adaptive thresholding for the maximum
                #   - estimate probability density of the maximum

                tr2 = np.mean(intensity_v)*0.25
                max_pks = find_peaks(sig+min(sig),prominence=tr2,distance=hr_win)[0]

        return max_pks

    ###########################################################################
    ############################ Bandpass filtering ###########################
    ###########################################################################
    def def_bandpass(self, sig: np.array, fs: int, lower_cutoff: float, upper_cutoff: float):
        """
        def_bandpass filter function detects all peaks in the raw and also in the filtered signal to find.

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array
        :param fs: sampling frequency
        :type fs: int
        :param lower_cutoff: lower cutoff frequency
        :type lower_cutoff: float
        :param upper_cutoff: upper cutoff frequency
        :type upper_cutoff: float

        :return: bandpass filtered signal, 1-d array.

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
    def elim_vlfs(self, s: np.array, up: DotMap, lower_cutoff: float):
        """
        This function filter the high frequency components.

        :param s: array of signal with shape (N,) where N is the length of the signal
        :type s: 1-d array
        :param up: setup up parameters of the algorithm
        :type up: DotMap
        :param lower_cutoff: lower cutoff frequency
        :type lower_cutoff: float

        :return: high frequency filtered signal, 1-d array.

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
        AMfilter = b

        s_filt=DotMap()
        try:
            s_filt.v = filtfilt(AMfilter, 1, s.v)
            s_filt.v = s.v-s_filt.v
        except:
            s_filt.v = s.v

        s_filt.fs = s.fs

        return s_filt

    ###########################################################################
    ################### Filter the low frequency components  ##################
    ###########################################################################
    def elim_vhfs(self, s: np.array, up: DotMap, upper_cutoff: float):
        """
        This function filter the high frequency components.

        :param s: array of signal with shape (N,) where N is the length of the signal
        :type s: 1-d array
        :param up: setup up parameters of the algorithm
        :type up: DotMap
        :param upper_cutoff: lower cutoff frequency
        :type upper_cutoff: float

        :return: low frequency filtered signal, 1-d array.

        """


        ## Filter signal to remove VHFs
        # Adapted from RRest
        s_filt = DotMap()

        ##Eliminate nans
        s.v[np.isnan(s.v)] = np.mean(s.v[~np.isnan(s.v)])

        ##Check to see if sampling freq is at least twice the freq of interest
        if (up.paramSet.elim_vhf.Fpass/(s.fs/2)) >= 1:
            s_filt.v = s.v
            return

        fc = upper_cutoff
        ripple = -20 * np.log10(up.paramSet.elim_vhf.Dstop)
        width = abs(up.paramSet.elim_vhf.Fpass - up.paramSet.elim_vhf.Fstop) / (s.fs / 2)
        [N, beta] = kaiserord(ripple, width)
        if N * 3 > len(s):
            N = round(N / 3)
        b = firwin(N, fc * 2 / s.fs, window=('kaiser', beta), scale=('True'))
        AMfilter = b

        ## Remove VHFs
        s_dt=detrend(s.v)
        s_filt.v = filtfilt(AMfilter, 1, s_dt)

        return s_filt

    ###########################################################################
    ########################### Heart Rate estimation #########################
    ###########################################################################
    def estimate_HR(self, sig: np.array, fs: int, up: DotMap, hr_past: int):
        """
        Heart Rate Estimation function estimate the heart rate according to the previous heart rate in given time window

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array
        :type fs: int
        :param up: setup up parameters of the algorithm
        :type up: DotMap
        :param hr_past: the average heart rate in the past in given time window
        :type hr_past: int

        :return: estimated heart rate, 1-d array.

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
        hr = int(hr[0])

        return hr

    ###########################################################################
    ############# Estimate derivative from highly filtered signal #############
    ###########################################################################
    def estimate_deriv(self, sig: np.array):
        """
        Derivative Estimation function estimate derivative from highly filtered signal based on the
        General least-squares smoothing and differentiation by the convolution (Savitzky Golay) method

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array

        :return: derivative, 1-d array.

        """

        #Savitzky Golay
        deriv_no = 1
        win_size = 5
        deriv = self.savitzky_golay(sig, deriv_no, win_size)

        return deriv


    def savitzky_golay(self, sig: np.array, deriv_no: int, win_size: int):
        """
        This function estimate the Savitzky Golay derivative from highly filtered signal

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array
        :param deriv_no: number of derivative
        :type deriv_no: int
        :param win_size: size of window
        :type win_size: int

        :return: Savitzky Golay derivative, 1-d array.

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
    def find_pulse_peaks(self, p2: np.array, p3: np.array):
        """
        Pulse detection function detect the pulse peaks according to the peaks of 1st and 2nd derivatives
        General least-squares smoothing and differentiation by the convolution (Savitzky Golay) method

        :param p2: peaks of the 1st derivatives
        :type p2: 1-d array
        :param p3: peaks of the 2nd derivatives
        :type p2: 1-d array

        :return: pulse peaks, 1-d array.

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
    def  correct_IBI(self, p: np.array, m: np.array, hr: float, fs: int, up: DotMap):
        """
        This function corrects the peaks' location (interbeat intervals) error

        :param p: systolic peaks of the PPG signal
        :type p: 1-d array
        :param m: all maxima of the PPG signal
        :type m: 1-d array
        :param hr: heart rate
        :type hr: float
        :param fs: sampling frequency
        :type fs: int
        :param up: setup up parameters of the algorithm
        :type up: DotMap

        :return: onsets, 1-d array.

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
        fp = self.find_reduced_IBIs(d, hr, up)

        # remove FPs
        pc = np.array(pc1)[fp]

        # Correct false negatives
        d = np.diff(pc)/fs    # interbeat intervals in secs
        fn = self.find_prolonged_IBIs(d, hr, up)

        pc = pc1

        return pc, fn

    def find_reduced_IBIs(self,IBIs: np.array, med_hr: float, up: DotMap):
        """
        This function finds the reduced interbeat intervals

        :param IBIs: interbeat intervals in secs
        :type IBIs: 1-d array
        :param med_hr: median heart rate
        :type med_hr: float
        :param up: setup up parameters of the algorithm
        :type up: DotMap

        :return: fp, the false positive case

        """

        IBI_thresh = up.lower_hr_thresh_prop*60/med_hr
        fp = IBIs < IBI_thresh
        fp = [*np.where(fp == 0)[0].astype(int)]
        return fp

    def find_prolonged_IBIs(self, IBIs: np.array, med_hr: float, up: DotMap):
        """
        This function finds the prolonged interbeat intervals

        :param IBIs: interbeat intervals in secs
        :type IBIs: 1-d array
        :param med_hr: median heart rate
        :type med_hr: float
        :param up: setup up parameters of the algorithm
        :type up: DotMap

        :return: fn, the false negative case

        """
        IBI_thresh = up.upper_hr_thresh_prop*60/med_hr
        fn = IBIs > IBI_thresh

        return fn

    ###########################################################################
    ####################### Setup up the beat detector ########################
    ###########################################################################
    def set_beat_detection(self):
        """
        This function setups the filter parameters of the algorithm

        :return: filter parameters of the algorithm, DotMap

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
    def find_onsets(self,sig: np.array, fs: int, up: DotMap, peaks: np.array, med_hr: float):
        """
        This function finds the onsets of PPG sigal

        :param sig: array of signal with shape (N,) where N is the length of the signal
        :type sig: 1-d array
        :param fs: sampling frequency
        :type fs: int
        :param up: setup up parameters of the algorithm
        :type up: DotMap
        :param peaks: peaks of the signal
        :type peaks: 1-d array
        :param med_hr: median heart rate
        :type med_hr: float

        :return: onsets, 1-d array

        """

        Y1=self.def_bandpass(sig, fs, 0.9*up.fl_hz, 3*up.fh_hz)
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
    def get_dicrotic_notch (self, peaks: np.array, onsets: list):
        """
        Dicrotic Notch function estimate the location of dicrotic notch in between the systolic and diastolic peak

        :param peaks: peaks of the signal
        :type peaks: 1-d array
        :param onsets: onsets of the signal
        :type onsets: list

        :return: location of dicrotic notches, 1-d array
        """

        ## The 2nd derivative and Hamming low pass filter is calculated.
        dxx = np.diff(np.diff(self.ppg))
        fs = self.fs

        # Make filter
        Fn = fs / 2                                 # Nyquist Frequency
        FcU = 20                                    # Cut off Frequency: 20 Hz
        FcD = FcU + 5                               # Transition Frequency: 5 Hz

        n = 21                                      # Filter order
        f = [0, (FcU / Fn), (FcD / Fn), 1]          # Frequency band edges
        a = [1, 1, 0, 0]                            # Amplitudes
        b = firls(n, f, a)

        lp_ppg = filtfilt(b, 1,  dxx)          # Low pass filtered signal with 20 cut off Frequency and 5 Hz Transition width

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
            nth_beat = lp_ppg[onsets[i]:onsets[i + 1]]

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
            ## NOTE!! Shifting with 26 ms. FIX IT!
            shift=int(self.fs*0.026)
            dic_not.append(max_pp_i+onsets[i]+shift)

        return dic_not

    ###########################################################################
    ########################## Detect diastolic peak ##########################
    ###########################################################################
    def get_diastolic_peak(self, onsets: list, dicroticnotch: list, e_point: pd.Series):
        """
        Dicrotic Notch function estimate the location of dicrotic notch in between the systolic and diastolic peak

        :param onsets: onsets of the signal
        :type onsets: list
        :param dicroticnotches: dicrotic notches of the signal
        :type dicroticnotches: list
        :param e_point: e-points of the signal
        :type e_point: pd.Series

        :return: diastolicpeak location of dicrotic notches, 1-d array
        """

        nan_v = np.empty(len(dicroticnotch))
        nan_v[:] = np.NaN
        diastolicpeak = nan_v

        for i in range(0,len(dicroticnotch)):
            try:
                len_segments=np.diff(onsets)*0.80
                end_segment=int(onsets[i]+len_segments[i])
                try:
                    start_segment = int(dicroticnotch[i])
                    temp_segment = self.ppg[start_segment:end_segment]
                    max_locs, _ = find_peaks(temp_segment)

                    if len(max_locs)==0:
                        start_segment = int(e_point[i])
                        temp_segment = self.vpg[start_segment:end_segment]
                        max_locs, _ = find_peaks(temp_segment)

                except:
                    pass

                max_dn = max_locs[0] + start_segment
                diastolicpeak[i] = max_dn
            except:
                pass

        return diastolicpeak

    ###########################################################################
    ####################### Get First Derivitive Points #######################
    ###########################################################################
    def get_vpg_fiducials(self, onsets: list):
        """Calculate first derivitive points u and v from the PPG' signal

        :param onsets: onsets of the signal
        :type onsets: list

        :return:
            - u: The highest amplitude between the pulse onset and systolic peak on PPG'
            - v: The lowest amplitude between the u-point and diastolic peak on PPG'
            - w: The first local maximum or inflection point after the dicrotic notch on PPG’
        """

        dx = self.vpg

        nan_v = np.empty(len(onsets)-1)
        nan_v[:] = np.NaN
        u, v, w = copy.deepcopy(nan_v),copy.deepcopy(nan_v),copy.deepcopy(nan_v),

        for i in range(0,len(onsets)-1):
            try:
                segment = dx[onsets[i]:onsets[i + 1]]

                # u fiducial point
                max_loc = np.argmax(segment)+onsets[i]
                u[i]=max_loc

                # v fiducial point
                upper_bound_coeff = 0.66
                v_upper_bound = ((onsets[i + 1] - onsets[i]) * upper_bound_coeff + onsets[i]).astype(int)
                min_loc = np.argmin(dx[int(u[i]):v_upper_bound])+u[i]-1
                v[i] = min_loc

                # w fiducial point
                temp_segment=dx[int(v[i]):onsets[i+1]]
                max_locs, _ = find_peaks(temp_segment)
                max_w = max_locs[0] + v[i] - 1
                w[i] = max_w

            except:
                pass

        vpg_fp = pd.DataFrame({"u":[], "v":[], "w":[]})
        vpg_fp.u, vpg_fp.v, vpg_fp.w = u, v, w
        return vpg_fp

    ###########################################################################
    ####################### Get Second Derivitive Points ######################
    ###########################################################################
    def get_apg_fiducials(self, onsets: list, peaks: np.array):
        """Calculate Second derivitive points a, b, c, d, e, and f from the PPG" signal

        :param onsets: onsets of the signal
        :type onsets: list
        :param peaks: peaks of the signal
        :param types: 1-d array

        :return:
            - a: The highest amplitude between pulse onset and systolic peak on PPG"
            - b: The first local minimum after the a-point on PPG"
            - c: The local maximum with the highest amplitude between the b-point and e-point, or if no local maximum is present then the inflection point on PPG"
            - d: The local minimum with the lowest amplitude between the c-point and e-point, or if no local minimum is present then the inflection point on PPG"
            - e: The local maximum with the highest amplitude after the b-point and before the diastolic peak on PPG"
            - f: The first local minimum after the e-point on PPG"
        """

        sig = self.ppg
        ddx = self.apg
        dddx = self.jpg

        nan_v = np.empty(len(onsets)-1)
        nan_v[:] = np.NaN
        a, b, c, d, e, f = copy.deepcopy(nan_v),copy.deepcopy(nan_v),copy.deepcopy(nan_v),copy.deepcopy(nan_v),copy.deepcopy(nan_v),copy.deepcopy(nan_v)
        for i in range(0,len(onsets)-1):

            try:
                # a fiducial point
                temp_pk=np.argmax(sig[onsets[i]:onsets[i + 1]])+onsets[i]-1
                temp_segment=ddx[onsets[i]:temp_pk]
                max_locs, _ = find_peaks(temp_segment)
                try:
                    max_loc = max_locs[np.argmax(temp_segment[max_locs])]
                except:
                    max_loc = temp_segment.argmax()

                max_a = max_loc + onsets[i] - 1
                a[i] = max_a

                # b fiducial point
                temp_segment=ddx[int(a[i]):onsets[i+1]]
                min_locs, _ = find_peaks(-temp_segment)
                min_b = min_locs[0] + a[i] - 1
                b[i] = min_b

                # e fiducial point
                e_lower_bound = peaks[i]
                upper_bound_coeff = 0.6
                e_upper_bound = ((onsets[i + 1] - onsets[i]) * upper_bound_coeff + onsets[i]).astype(int)
                temp_segment=ddx[int(e_lower_bound):int(e_upper_bound)]
                max_locs, _ = find_peaks(temp_segment)
                if max_locs.size==0:
                    temp_segment=ddx[int(peaks[i]):onsets[i + 1]]
                    max_locs, _ = find_peaks(temp_segment)

                try:
                    max_loc = max_locs[np.argmax(temp_segment[max_locs])]
                except:
                    max_loc = temp_segment.argmax()

                max_e = max_loc + e_lower_bound - 1
                e[i] = max_e

                # c fiducial point
                temp_segment = ddx[int(b[i]):int(e[i])]
                max_locs, _ = find_peaks(temp_segment)
                if max_locs.size>0:
                    max_loc = max_locs[0]
                else:
                    temp_segment = dddx[int(b[i]):int(e[i])]
                    min_locs, _ = find_peaks(-temp_segment)

                    if min_locs.size > 0:
                        max_loc = min_locs[np.argmin(temp_segment[min_locs])]
                    else:
                        max_locs, _ = find_peaks(temp_segment)
                        try:
                            max_loc = max_locs[0]
                        except:
                            max_loc = temp_segment.argmax()

                max_c = max_loc + b[i] - 1
                c[i] = max_c

                # d fiducial point
                temp_segment = ddx[int(c[i]):int(e[i])]
                min_locs, _ = find_peaks(-temp_segment)
                if min_locs.size > 0:
                    min_loc = min_locs[np.argmin(temp_segment[min_locs])]
                    min_d = min_loc + c[i] - 1
                else:
                    min_d = max_c

                d[i] = min_d

                # f fiducial point
                temp_segment = ddx[int(e[i]):onsets[i + 1]]
                min_locs, _ = find_peaks(-temp_segment)
                if (min_locs.size > 0) and (min_locs[0]<len(sig)*0.8):
                    min_loc = min_locs[0]
                else:
                    min_loc = 0

                min_f = min_loc + e[i] - 1
                f[i] = min_f
            except:
                pass

        apg_fp = pd.DataFrame({"a":[], "b":[], "c":[],"d":[], "e":[], "f":[]})
        apg_fp.a, apg_fp.b, apg_fp.c, apg_fp.d, apg_fp.e, apg_fp.f = a, b, c, d, e, f
        return apg_fp

    def get_jpg_fiducials(self, onsets: list, apg_fp: pd.DataFrame):
        """Calculate third derivitive points p1 and p2 from the PPG'" signal

            :param onsets: onsets of the signal
            :type onsets: list
            :param apg_fp: fiducial points of PPG" signal
            :type apg_fp: DataFrame

            :return:
                - p1: The first local maximum after the b-point on PPG'"
                - p2: The last local minimum after the b-point and before the d-point on PPG'"

        """
        dddx = self.jpg

        nan_v = np.empty(len(onsets)-1)
        nan_v[:] = np.NaN
        p1, p2 = copy.deepcopy(nan_v), copy.deepcopy(nan_v)

        for i in range(0, len(onsets) - 1):
            try:
                # p1 fiducial point
                ref_b = apg_fp.b[np.squeeze(np.where(np.logical_and(apg_fp.b > onsets[i], apg_fp.b < onsets[i + 1])))]
                if ref_b.size == 0:
                    ref_b = onsets[i]

                temp_segment = dddx[int(ref_b):onsets[i + 1]]
                max_locs, _ = find_peaks(temp_segment)
                try:
                    max_loc = max_locs[0]
                except:
                    max_loc = temp_segment.argmax()

                max_p1 = max_loc + ref_b - 1
                p1[i] = max_p1

                # p2 fiducial point
                ref_start = p1[i]
                ref_c = apg_fp.c[np.squeeze(np.where(np.logical_and(apg_fp.c > onsets[i], apg_fp.c < onsets[i + 1])))]
                ref_d = apg_fp.d[np.squeeze(np.where(np.logical_and(apg_fp.d > onsets[i], apg_fp.d < onsets[i + 1])))]

                if ref_d > ref_c:
                    ref_end = ref_d
                    min_ind = -1
                elif ref_c.size > 0:
                    ref_start = ref_c
                    ref_end = onsets[i + 1]
                    min_ind = 0
                elif apg_fp.e.size > 0:
                    ref_end = onsets[i + 1]
                    min_ind = 0

                temp_segment = dddx[int(ref_start):int(ref_end)]
                min_locs, _ = find_peaks(-temp_segment)
                if min_locs.size > 0:
                    min_p2 = min_locs[min_ind] + ref_start - 1
                else:
                    min_p2 = ref_end

                p2[i] = min_p2

            except:
                pass

        jpg_fp = pd.DataFrame({"p1": [], "p2": []})
        jpg_fp.p1, jpg_fp.p2 = p1, p2

        return jpg_fp

    def correct_fiducials(self,fiducials=pd.DataFrame(), correction=pd.DataFrame()):
        """Correct the Fiducial Points

            :param fiducials: DataFrame where the key is the name of the fiducial points and the value is the list of fiducial points PPG Fiducials Points
            :type fiducials: DataFrame
            :param correction: DataFrame where the key is the name of the fiducial points and the value is bool
            :type correction: DataFrame

            :return:
                - fiducials: a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points
        """

        for i in range(0,len(fiducials.on)):

            # Correct onset
            if correction.on[0]:
                try:
                    win_onr = self.fs * 0.005
                    if fiducials.on[i]>win_onr:
                        win_onl = win_onr
                    else:
                        win_onl = fiducials.on[i]

                    min_loc = np.argmax(self.ppg[fiducials.on[i]+win_onl:fiducials.on[i]+win_onr]) + fiducials.on[i]
                    if fiducials.on[i] != min_loc:

                        if fiducials.a[i] > self.fs*0.075:
                            win_a = int(self.fs*0.075)
                        else:
                            win_a = int(fiducials.a[i])

                        fiducials.on[i] = np.argmax(self.jpg[int(fiducials.a[i]) - win_a:int(fiducials.a[i])]) + int(fiducials.a[i]) - win_a
                except:
                    pass

            # Correct dicrotic notch
            if correction.dn[0]:
                try:
                    temp_segment = self.ppg[int(fiducials.sp[i]):int(fiducials.dp[i])]
                    min_dn = find_peaks(-temp_segment)[0] + fiducials.sp[i]
                    diff_dn = abs(min_dn - fiducials.dp[i])
                    if len(min_dn) > 0 and diff_dn > round(self.fs / 100):
#                        fiducials.dn[i] = min_dn
                        fiducials.loc[i, 'dn'] = min_dn
                        try:
                            strt_dn = int(fiducials.sp[i])
                            stp_dn = int(fiducials.f[i])
#                            fiducials.dn[i] = find_peaks(-self.ppg[strt_dn:stp_dn])[0][-1] + strt_dn
                            fiducials.loc[i, 'dn'] = find_peaks(-self.ppg[strt_dn:stp_dn])[0][-1] + strt_dn
                            if fiducials.dn[i] > min_dn:
#                                fiducials.dn[i] = min_dn
                                fiducials.loc[i, 'dn'] = min_dn
                        except:
                            strt_dn = fiducials.e[i]
                            stp_dn = fiducials.f[i]
                            fiducials.dn[i] = np.argmin(self.jpg[strt_dn:stp_dn]) + strt_dn
                            if fiducials.dn[i] > min_dn:
#                                fiducials.dn[i] = min_dn
                                fiducials.loc[i, 'dn'] = min_dn

                except:
                    pass

            # Correct w-point
            if correction.w[0]:
                if fiducials.w[i] > fiducials.f[i]:
                    fiducials.loc[i, 'w'] = fiducials.f[i]

                if fiducials.w[i] < fiducials.e[i]:
                    try:
                        fiducials.loc[i, 'w'] = [np.argmax(self.vpg[int(fiducials.e[i]):int(fiducials.f[i])]) + fiducials.e[i]]
                    except:
                        pass

            # Correct v-point and w-point
            if correction.v[0] and correction.w[0]:
                if fiducials.v[i] > fiducials.e[i]:
                    try:
                        fiducials.loc[i, 'v'] = [np.argmin(self.vpg[int(fiducials.u[i]):int(fiducials.e[i])]) + fiducials.u[i]]
                        fiducials.loc[i, 'w'] = [find_peaks(self.vpg[int(fiducials.v[i]):int(fiducials.f[i])])[0][0] + fiducials.v[i]]
                    except:
                        pass

            # Correct f-point
            if correction.f[0]:
                try:
                    temp_end=int(np.diff(fiducials.on[i:i+2])*0.8)
                    temp_segment=self.apg[int(fiducials.e[i]):int(fiducials.on[i]+temp_end)]
                    min_f=np.argmin(temp_segment)+fiducials.e[i]

                    if fiducials.w[i] > fiducials.f[i]:
                        fiducials.f[i] = min_f
                except:
                    pass

        # Correct diastolic peak
        if correction.dp[0]:
            try:
                fiducials.dp = self.get_diastolic_peak(fiducials.on, fiducials.dn, fiducials.e)
            except:
                pass

        return fiducials
