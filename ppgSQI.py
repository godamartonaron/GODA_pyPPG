import numpy as np
from scipy.signal import detrend, find_peaks, correlate
from numpy.fft import fft, ifft, fftshift

def ppgSQI(ppg, fs, ann_ppg):
    # PPG Signal Quality Index based on beat template correlation.
    # input:
    #   ppg: PPG data
    #   fs: Samples frequency
    #   annppg: PPG annotation time(samples)
    # output:
    #   psqi: PPG SQI
    #
    # Author:
    # Marton A. Goda: Faculty of Biomedical Engineering,
    # Technion – Israel Institute of Technology, Haifa, Israel(October 2022)
    #
    # Original Matlab implementation:
    # Qiao Li, November 10, 2014.
    # https://github.com/MIT-LCP/PhysioNetChallengePublic/blob/master/2015/sample-submission/ppgSQI.m

    Fs = fs
    # Create PPG template
    t,v = template(ppg, ann_ppg-1, fs)

    c1 = np.empty(len(ann_ppg)-1)
    c1[:] = np.NaN
    psqi = np.empty(len(ann_ppg)-1)
    psqi[:] = np.NaN

    for j in range (0,len(ann_ppg) - 1):
        # Calculate correlation coefficients based on the template length
        beatbegin = ann_ppg[j]-1
        beatend = ann_ppg[j + 1]-1
        if beatend - beatbegin > 3 * Fs:
            beatend = beatbegin + 3 * Fs

        templatelength = len(t)
        if (beatbegin + templatelength - 1 > len(ppg)) or (beatend > len(ppg)) or (beatbegin < 1):
            continue

        currentb = j
        cc = np.corrcoef(t, ppg[beatbegin:beatbegin + templatelength])
        c1[j] = cc[0, 1]
        if (c1[j] < 0):
            c1[j] = 0

        psqi[currentb] = c1[j]

    return psqi

def template(wave, anntime, samp_freq):
    # PPG waveform template creation.
    # Written by Qiao Li, February 21, 2011.
    #
    # input:
    #   wave: PPG data
    #   anntime: PPG annotation time(sample)
    #   samp_freq: sampling frequency, default is 125 Hz
    # output:
    #   t: PPG waveform template based on normal - length beats
    #   valid:  1 for valid template
    #           0 for invalid template

    t = []
    valid = 0

    # according to heart rate max(300 bpm) and min(20 bpm) to get max and min beat - by - beat interval
    hr_max = 300
    bb_interval_min = samp_freq * 60 / hr_max
    hr_min = 20
    bb_interval_max = samp_freq * 60 / hr_min

    # Normal beat thresholds
    normal_beat_length_min = 0.7
    normal_beat_lentth_max = 1.5
    normal_beat_percent_threshold = 0.5

    # using correlate to get the basic period of the PPG as the length of template
    data = detrend(wave)

    y = correlate(data, data, 'full', method='fft')
    lenw = len(wave)
    lena = len(anntime)
    i = lenw

    locs = find_peaks(y[i:])[0]
    pks = y[locs]
    if len(pks)==0:
        return

    max_v=max(pks)
    max_i=np.where(pks==max_v)[0][0]
    i = locs[max_i]

    cycle = samp_freq
    if i < lenw - 1:
        cycle = i

    # cumulate the beats with reasonable length to get template
    if lena < 2:
        return

    p0 = 1
    i = anntime[p0]

    temp_ahead=0
    while i - temp_ahead < 1:
        p0 = p0 + 1
        if (p0 > lena):
            t = wave
            valid = 0
            return

        i = anntime[p0]

    if p0 + 1 >= lena:
        return

    beat_interval = np.diff(anntime[p0:len(anntime)])
    median_bi = np.median(beat_interval)

    if median_bi!='NaN':
        temp_peak = abs(locs - median_bi)
        m = min(temp_peak)
        i = np.where(temp_peak==m)[0][0]
        cycle = locs[i]+1
    else:
        return

    # the length of template valid detection
    valid = 1
    if (cycle > bb_interval_max) or (cycle < bb_interval_min):
        valid = 0
        t = np.zeros(cycle)
        return

    n = 0
    d1 = 0
    invalidn = 0
    currentbeatlength = anntime[p0 + 1] - anntime[p0]
    if currentbeatlength > 0:
        d1 = wave[i - temp_ahead:i + cycle]
        n = 1
    else:
        invalidn = invalidn + 1
        d1 = np.zeros(cycle + temp_ahead)

    p0 = p0
    if p0 < lena - 1:
        i = anntime[p0]
        n = 1
        invalidn = 0
        while (i < lenw - cycle) and (p0 < lena - 1):
            currentbeatlength = anntime[p0 + 1] - anntime[p0]
            if currentbeatlength > 0:
                d1 = d1 + wave[i - temp_ahead:i + cycle]
                n = n + 1
            else:
                invalidn = invalidn + 1

            p0 = p0 + 1
            i = anntime[p0]

        d1 = d1/n
        # normal beat is less than the reasonable percentage of all beats
        if (n / (n + invalidn)) < normal_beat_percent_threshold:
            valid = 0
        else:
            valid = 0
    else:
        valid=0

    t = d1

    return t,valid