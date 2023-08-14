import pyPPG
from pyPPG.preproc import Preprocessing

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import os

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
def load_data(data_path: str, start = 0, end = 0, filtering = True):
    """
    Load PPG data function load the raw PPG data.

    :param data_path: path of the PPG signal
    :type data_path: str
    :param start: beginning the of signal in sample
    :type start: int
    :param end: end of the signal in sample
    :type end: int
    :param filtering: a bool for filtering
    :type filtering: bool

    :return s: a struct of PPG signal:

        * s.start: beginning of the signal in sample
        * s.end: end of the signal in sample
        * s.v: a vector of PPG values
        * s.fs: the sampling frequency of the PPG in Hz
        * s.name: name of the record
        * s.v: 1-d array, a vector of PPG values
        * s.fs: the sampling frequency of the PPG in Hz
        * s.filt_sig: 1-d array, a vector of the filtered PPG values
        * s.filt_d1: 1-d array, a vector of the filtered PPG' values
        * s.filt_d2: 1-d array, a vector of the filtered PPG" values
        * s.filt_d3: 1-d array, a vector of the filtered PPG'" values
    """

    if data_path=="":
        sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl .txt")])
    else:
        sig_path=data_path

    if sig_path.rfind('/')>0:
        start_c = sig_path.rfind('/')+1
    else:
        start_c = sig_path.rfind('\\')

    stop_c=sig_path.rfind('.')
    rec_name=sig_path[start_c:stop_c]

    try:
        sig_format=sig_path[len(sig_path)-sig_path[::-1].index('.'):]
    except:
        print('Invalid signal path!')

    if sig_format=='mat':
        input_sig = scipy.io.loadmat(sig_path)
        sig = np.float64(np.squeeze(input_sig.get("Data")))[0:]
        try:
            fs = np.squeeze(input_sig.get("Fs"))
        except:
            fs = 100
            print('The default sampling frequency is 100 Hz for .mat.')
    elif sig_format=='csv':
        input_sig = np.loadtxt(sig_path, delimiter=',').astype(int)
        sig = input_sig
        fs = 75
        print('The default sampling frequency is 75 Hz for .csv.')
    elif sig_format=='txt':
        try:
            input_sig = np.loadtxt(sig_path, delimiter='\t').astype(int)
        except:
            try:
                input_sig = np.loadtxt(sig_path, delimiter=' ').astype(int)
            except:
                print('ERROR! The data separator is not supported for .txt.')
        sig = input_sig
        fs = 1000
        print('The default sampling frequency is 1 kHz for .txt.')
    elif sig_format == 'edf':
        input_sig = mne.io.read_raw_edf(sig_path)
        sig = -input_sig['Pleth'][0][0]
        try:
            fs = int(np.round(input_sig.info['sfreq']))
        except:
            fs = 256
            print('The default sampling frequency is 256 Hz for .edf.')

    s = DotMap()

    s.start = start
    if start<end:
        s.end = end
    else:
        s.end = len(sig)

    s.v=sig[s.start:s.end]
    s.fs=fs
    s.name=rec_name
    s.filt_sig, s.filt_d1, s.filt_d2, s.filt_d3 = Preprocessing(s, filtering = True)

    return s

###########################################################################
########################### Plot Fiducial points ##########################
###########################################################################
def plot_fiducials(s: pyPPG.PPG, fp: pyPPG.Fiducials, savingfolder: str):
    """
    Plot fiducial points of the filtered PPG signal.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :param fp: a struct of fiducial points
    :type fp: pyPPG.Fiducials object
    :param savingfolder: location of the saved figure
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

    major_ticks_names = range(0, int(len_sig/s.fs),step_big)
    len_ticks_names=len(major_ticks_names)
    major_diff=len_sig/len_ticks_names
    minor_diff = len_sig / len_ticks_names / step_big
    major_ticks = np.arange(str_sig, end_sig, major_diff)
    minor_ticks = np.arange(str_sig, end_sig, minor_diff)

    for n in fid_names:
        ind = fid_names.index(n)
        if s_type[ind][-1] == 'g':
            plt_num = 1
            plt.subplot(411)
            plt.title(s.name, fontsize=20)
        else:
            plt_num = int(s_type[ind][-1]) + 1


        ax = plt.subplot(4, 1, plt_num)
        tmp_pnt=eval("fp." + n + ".values")
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
    plt.xticks(major_ticks,major_ticks_names, fontsize=20)
    plt.show()

    canvas = FigureCanvas(fig)
    tmp_dir=savingfolder+os.sep+'PPG_Figures'+os.sep

    os.makedirs(tmp_dir, exist_ok=True)

    canvas.print_png((tmp_dir+'%s_btwn_%s-%s.png') % (s.name,s.start,s.end))
    print('Figure has been saved in the "'+savingfolder+'".')

###########################################################################
################################# Save Data ###############################
###########################################################################

def save_data(s: pyPPG.PPG, fp: pyPPG.Fiducials, bm: pyPPG.Biomarkers, savingformat: str, savingfolder: str):
    """
    Save the results of the filtered PPG analysis.

    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :type s: DotMap
    :param fp: a struct of fiducial points
    :type fp: pyPPG.Fiducial object
    :param bm: a dictionary of biomarkers
    :type bm: pyPPG.Biomarkers object
    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str
    :param savingfolder: location of the saved data
    :type savingfolder: str
    """

    tmp_dir = savingfolder
    os.makedirs(tmp_dir, exist_ok=True)

    temp_dirs = ['Fiducial_points', 'Biomarker_vals', 'Biomarker_stats', 'Biomarker_defs', 'PPG_struct']
    for i in temp_dirs:
        temp_dir = tmp_dir + os.sep + i + os.sep
        os.makedirs(temp_dir, exist_ok=True)

    keys=s.__dict__.keys()
    keys_list = list(keys)
    sc=DotMap()
    for i in keys_list:
        exec('sc.'+i+' = s.'+i)

    file_name = (r'.' + os.sep + tmp_dir + os.sep + temp_dirs[4] + os.sep + s.name + '_data_btwn_%s-%s.mat')%(s.start,s.end)
    scipy.io.savemat(file_name, sc)

    BM_keys = bm.bm_vals.keys()

    if savingformat=="csv":
        file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[0]+os.sep+s.name+'_'+'Fiducials_btwn_%s-%s.csv')%(s.start,s.end)
        fp.get_fp().to_csv(file_name)

        for key in BM_keys:
            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[1]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start,s.end)
            bm.bm_vals[key].to_csv(file_name,index=True,header=True)

            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[2]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start,s.end)
            bm.bm_stats[key].to_csv(file_name, index=True, header=True)

            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[3]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start,s.end)
            bm.bm_defs[key].to_csv(file_name, index=True, header=True)

    elif savingformat=="mat":
        matlab_struct = fp.get_fp().to_dict(orient='list')
        file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[0]+os.sep+s.name+'_'+'Fiducials_btwn_%s-%s.mat')%(s.start,s.end)
        scipy.io.savemat(file_name,matlab_struct)

        for key in BM_keys:

            matlab_struct = bm.bm_vals[key].to_dict(orient='list')
            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[1]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start,s.end)
            scipy.io.savemat(file_name,matlab_struct)

            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[2]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start,s.end)
            scipy.io.savemat(file_name, bm.bm_stats[key])

            matlab_struct = bm.bm_defs[key].to_dict(orient='list')
            file_name = (r'.'+os.sep+tmp_dir+os.sep+temp_dirs[3]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start,s.end)
            scipy.io.savemat(file_name,matlab_struct)
    else:
        print('The file format is not suported for data saving! You can use "mat" or "csv" file formats.')

    print('Results have been saved into the "'+tmp_dir+'".')