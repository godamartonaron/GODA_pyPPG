from Prefiltering import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time


def load_data(filtering):
    # sig_path = 'D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/GODA_pyPPG/sample_data/PPG_sample_00.mat'
    sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl")])

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

    if filtering:
        s.filt_sig, s.filt_d1, s.filt_d2, s.filt_d3 = Prefiltering(s)

    return s


def plot_fiducials(s,fiducials):
    fig = plt.figure(figsize=(8, 9))
    ax1 = plt.subplot(411)
    plt.plot(s.filt_sig, 'k', label=None)
    ax2 = plt.subplot(412, sharex=ax1)
    plt.plot(s.filt_d1, 'k', label=None)
    ax3 = plt.subplot(413, sharex=ax2)
    plt.plot(s.filt_d2, 'k', label=None)
    ax4 = plt.subplot(414, sharex=ax3)
    plt.plot(s.filt_d3, 'k', label=None)
    fig.subplots_adjust(hspace=0, wspace=0)

    marker = ['o', 's', 's','o', 'o', 's', 'o', 'o', 's', 'o', 's', 'o', 's', 'o', 's']
    color = ['r', 'b', 'g','m', 'r', 'b', 'g', 'r', 'b', 'g', 'm', 'c', 'k', 'r', 'b']

    fid_names = ('sp', 'on', 'dn','dp', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')
    ylabe_names = ['PPG', 'PPG''', 'PPG\'\'', 'PPG\'\'\'']
    s_type = ['filt_sig', 'filt_sig', 'filt_sig','filt_sig', 'filt_d1', 'filt_d1', 'filt_d1', 'filt_d2', 'filt_d2', 'filt_d2', 'filt_d2', 'filt_d2', 'filt_d2', 'filt_d3','filt_d3']

    str_sig = 0
    end_sig = len(s.filt_sig)
    len_sig=end_sig-str_sig
    step_small = 1
    step_big = step_small * 5
    major_ticks = np.arange(str_sig, end_sig, int(step_big * s.fs))
    minor_ticks = np.arange(str_sig, end_sig, int(step_small * s.fs))

    for n in fid_names:
        ind = fid_names.index(n)
        if s_type[ind][-1] == 'g':
            plt_num = 1
            plt.subplot(411)
        else:
            plt_num = int(s_type[ind][-1]) + 1


        ax = plt.subplot(4, 1, plt_num)
        tmp_pnt=eval("fiducials['" + n + "'].values")
        tmp_pnt=tmp_pnt[~np.isnan(tmp_pnt)].astype(int)
        tmp_sig=eval("s." + s_type[ind])
        exec("plt.scatter(tmp_pnt, tmp_sig[tmp_pnt], s=60,linewidth=2, marker = marker[ind],facecolors='none', color=color[ind], label=n)")
        plt.ylabel(ylabe_names[plt_num - 1], fontsize=20)

        exec("plt.xlim([str_sig,end_sig])")
        leg = plt.legend(loc='upper right', fontsize=20, ncol=2,facecolor="orange")
        for text in leg.get_texts():
            text.set_weight('bold')

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)

        ax.grid(which='both')
        ax.grid(which='minor', alpha=0.2,axis='x')
        ax.grid(which='major', alpha=1,axis='x')

        plt.yticks([])


    plt.xlabel('Time [s]', fontsize=20)
    major_ticks_names = range(0, int(len_sig/s.fs),step_big)

    plt.xticks(major_ticks,major_ticks_names, fontsize=20)
    plt.show()

def save_data(fiducials,ppg_biomarkers,ppg_summary,ppg_statistics):
    fiducials.to_csv((r'./temp_dir/fidu/fiducial.csv'))
