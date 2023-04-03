from FiducialPoints import*
from Biomarkers2 import*
from Summary import*
from Statistics import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time

from six.moves import cPickle as pickle

import matplotlib.mlab

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':
    sig_path = 'D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/GODA_pyPPG/sample_data/PPG_sample_00.mat'
    #sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl")])

    sig_format=sig_path[len(sig_path)-sig_path[::-1].index('.'):]
    if sig_format=='mat':
        input_sig = scipy.io.loadmat(sig_path)
        hr = np.float64(np.squeeze(input_sig.get("Data")))[0:]
        fs = np.squeeze(input_sig.get("Fs"))
    elif sig_format=='csv':
        input_sig = np.loadtxt(sig_path, delimiter=',').astype(int)
        hr = input_sig
        fs = 75
    elif sig_format == 'edf':
        input_sig = mne.io.read_raw_edf(sig_path)
        hr=-input_sig[22][0][0]
        fs = 256

    s = DotMap()
    s.v=hr
    s.fs=fs

    Fid_time = time.time()
    fiducials = getFiducialsPoints(s.v,s.fs)

    #######
    ppg_biomarkers = Biomarkers2(s.v, s.fs, fiducials)
    ppg_summary = Summary(s.v, fiducials['peaks'], fiducials['onsets'], s.fs)
    ppg_statistics = Statistics(fiducials['peaks'], fiducials['onsets'], ppg_biomarkers)

    ##
    print('Detection time of fiducial points: '+str(round(time.time()-Fid_time,2))+' sec')

    plt.plot(hr,'k',linewidth=0.7)
    plt.plot(fiducials['peaks'], hr[fiducials['peaks']], 'ro')
    plt.plot(fiducials['onsets'], hr[fiducials['onsets']], 'bs')
    plt.plot(fiducials['dicroticnotch'], hr[fiducials['dicroticnotch']], 'm*')

    plt.legend(['sigal', 'peak','onset','dic.notch'])
    plt.xlabel('Sample (Fs='+str(s.fs)+' Hz)')
    plt.ylabel('Amplitude')
    plt.show()

    print('Program finished')
