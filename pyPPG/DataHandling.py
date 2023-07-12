from Prefiltering import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

###########################################################################
############################### Load PPG data #############################
###########################################################################
def load_data(filtering):
    """
    Load PPG data function load the raw PPG data.
    :param filtering: a bool for filtering
    :return: s: a struct of PPG signal:
            - s.v: a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.name: name of the record
            - s.filt_sig: a vector of PPG values
            - s.filt_d1: a vector of PPG values
            - s.filt_d2: a vector of PPG values
            - s.filt_d3: a vector of PPG values
    """

    sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl")])

    if sig_path.rfind('/')>0:
        start_c = sig_path.rfind('/')+1
    else:
        start_c = sig_path.rfind('\\')

    stop_c=sig_path.rfind('.')
    rec_name=sig_path[start_c:stop_c]


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
    s.name=rec_name

    if filtering:
        s.filt_sig, s.filt_d1, s.filt_d2, s.filt_d3 = Prefiltering(s)

    return s

###########################################################################
########################### Plot Fiducial points ##########################
###########################################################################
def plot_fiducials(s, fiducials, savefig):
    """
    Plot fiducial points of the filtered PPG signal.
    :param s: a struct of PPG signal:
            - s.v: a vector of PPG values
            - s.fs: the sampling frequency of the PPG in Hz
            - s.name: name of the record
            - s.filt_sig: a vector of PPG values
            - s.filt_d1: a vector of PPG values
            - s.filt_d2: a vector of PPG values
            - s.filt_d3: a vector of PPG values
    :param fiducials: a dictionary where the key is the name of the fiducial pints
            and the value is the list of fiducial points.
    :param savefig: a bool for fiducial points saving
    """

    fig = plt.figure(figsize=(20, 12))
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
    ylabe_names = ['PPG', 'PPG\'', 'PPG\'\'', 'PPG\'\'\'']
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
            plt.title(s.name, fontsize=20)
        else:
            plt_num = int(s_type[ind][-1]) + 1


        ax = plt.subplot(4, 1, plt_num)
        tmp_pnt=eval("fiducials['" + n + "'].values")
        tmp_pnt=tmp_pnt[~np.isnan(tmp_pnt)].astype(int)
        tmp_sig=eval("s." + s_type[ind])
        exec("plt.scatter(tmp_pnt, tmp_sig[tmp_pnt], s=60,linewidth=2, marker = marker[ind],facecolors='none', color=color[ind], label=n)")
        plt.ylabel(ylabe_names[plt_num - 1], fontsize=20)

        exec("plt.xlim([str_sig,end_sig])")
        leg = plt.legend(loc='upper right', fontsize=20, ncol=2,facecolor="orange",frameon=True)
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

    if savefig:
        canvas = FigureCanvas(fig)
        canvas.print_png(('temp_dir/PPG_Figures/%s.png') % (s.name))

###########################################################################
################################# Save Data ###############################
###########################################################################

def save_data(fiducials,ppg_biomarkers,ppg_statistics):
    """
    Save the results of the filtered PPG analysis.
    :param fiducials: a dictionary where the key is the name of the fiducial pints
            and the value is the list of fiducial points.
    :param biomarkers: dictionary of biomarkers in different categories:
        - PPG signal
        - Signal ratios
        - PPG derivatives
        - Derivatives ratios
    :param Statistics: data frame with summary of PPG features
    """

    fiducials.to_csv((r'./temp_dir/PPG_Fiducials/Fiducial.csv'))
    BM_keys=ppg_biomarkers.keys()
    for key in BM_keys:
        ppg_biomarkers[key].to_csv((r'./temp_dir/PPG_Biomarkers/%s.csv')% (key),index=True,header=True)
        ppg_statistics[key].to_csv((r'./temp_dir/PPG_Statistics/%s.csv') % (key), index=True, header=True)


