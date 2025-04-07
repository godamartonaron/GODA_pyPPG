import pyPPG

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import pandas as pd
from dotmap import DotMap
from tkinter import filedialog
import mne
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import os
import tkinter as tk
from tkinter import simpledialog
from scipy.io import savemat

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
def load_data(data_path = "", fs = np.nan, start_sig = 0, end_sig = -1, channel='Pleth', use_tk=True, print_flag=True):
    """
    Load raw PPG data.

    :param data_path: path of the file containing the PPG signal
    :type data_path: str
    :param start_sig: the first sample of the signal to be analysed
    :type start_sig: int
    :param fs: the sampling frequency of the PPG in Hz
    :type fs: int
    :param end_sig: the last sample of the signal to be analysed
    :type end_sig: int
    :param channel: channel of the .edf file
    :type channel: channel of the .edf file
    :param use_tk: a bool for using tkinter interface
    :type use_tk: bool
    :param print_flag: a bool for print message
    :type print_flag: bool

    :return: s: dictionary of the PPG signal:

        * s.start_sig: the first sample of the signal to be analysed
        * s.end_sig: the last sample of the signal to be analysed
        * s.v: a vector of PPG values
        * s.fs: the sampling frequency of the PPG in Hz
        * s.name: name of the record
        * s.v: 1-d array, a vector of PPG values
        * s.fs: the sampling frequency of the PPG in Hz
        * s.ppg: 1-d array, a vector of the filtered PPG values
        * s.vpg: 1-d array, a vector of the filtered PPG' values
        * s.apg: 1-d array, a vector of the filtered PPG" values
        * s.jpg: 1-d array, a vector of the filtered PPG'" values
        * s.filtering: a bool for filtering
        * s.correct: a bool for correcting fiducial points
    """

    if data_path=="":
        sig_path = filedialog.askopenfilename(title='Select SIGNAL file', filetypes=[("Input Files", ".mat .csv .edf .pkl .txt")])
    else:
        sig_path=data_path

    if sig_path.rfind('/')>0:
        start_c = sig_path.rfind('/')+1
    else:
        start_c = sig_path.rfind('\\')+1

    stop_c=sig_path.rfind('.')
    rec_name=sig_path[start_c:stop_c]

    try:
        sig_format=sig_path[len(sig_path)-sig_path[::-1].index('.'):]
    except:
        raise('Invalid signal path!')

    if sig_format=='mat':
        input_sig = scipy.io.loadmat(sig_path)
        sig = np.float64(np.squeeze(input_sig.get("Data")))[0:]
        try:
            fs = np.squeeze(input_sig.get("Fs"))
        except:
            fs = 125
            if print_flag: print('The default sampling frequency is 125 Hz for .mat.')
    elif sig_format=='csv':
        input_sig = pd.read_csv(sig_path, encoding='utf-8')
        sig = input_sig
        sig = np.squeeze(sig.values)

        if fs<=0:
            fs = 125
            if print_flag: print('The default sampling frequency is 125 Hz for .csv.')

    elif sig_format=='txt':
        try:
            input_sig = np.loadtxt(sig_path, delimiter='\t')
        except:
            try:
                input_sig = np.loadtxt(sig_path, delimiter=' ')
            except:
                print('ERROR! The data separator is not supported for .txt.')
        sig = input_sig

        if fs<=0:
            fs = 125
            if print_flag: print('The default sampling frequency is 125 Hz for .txt.')

    elif sig_format == 'edf':
        try:
            input_sig = mne.io.read_raw_edf(sig_path, include=channel)
            sig = -input_sig.get_data().squeeze()
        except:
            try:
                if use_tk: input_name = simpledialog.askstring("Input", "Define the PPG channel name:")
                input_sig = mne.io.read_raw_edf(sig_path, include=input_name)
                sig = -input_sig.get_data().squeeze()
            except:
                raise('There is no valid channel for PPG in the .edf file!')

        if fs<=0:
            try:
                fs = int(np.round(input_sig.info['sfreq']))
            except:
                fs = 256
                if print_flag: print('The default sampling frequency is 256 Hz for .edf.')

    s = DotMap()

    s.start_sig = start_sig
    if start_sig<end_sig:
        s.end_sig = end_sig
    else:
        s.end_sig = len(sig)

    try:
        s.v=sig[s.start_sig:s.end_sig]
    except:
        raise('There is no valid PPG signal!')

    s.fs=fs
    s.name=rec_name

    return s

###########################################################################
########################### Plot Fiducial points ##########################
###########################################################################
def plot_fiducials(s: pyPPG.PPG, fp: pyPPG.Fiducials, savefig=True, savingfolder='temp_dir', show_fig = True,
                   print_flag=True, use_tk=False, new_fig=True, marker=[], title='Detection', legend_loc='upper right',
                   legend_fontsize=20, marker_size=60, facecolor=False, subtext={}, canvas=np.nan):
    """
    Plot fiducial points of the filtered PPG signal.

    :param s: object of PPG signal
    :type s: pyPPG.PPG object
    :param fp: object of fiducial points
    :type fp: pyPPG.Fiducials object
    :param savefig: a bool for save figure
    :type savefig: bool
    :param savingfolder: location of the saved figure
    :param show_fig: a bool for show figure
    :type show_fig: bool
    :param print_flag: a bool for print message
    :type print_flag: bool
    :param use_tk: a bool for using tkinter interface
    :type use_tk: bool
    :param new_fig: a bool for creating new figure
    :type new_fig: bool
    :param marker: list of fiducial points markers
    :type marker: list
    :param title: title of the legend
    :type title: str
    :param legend_loc: location of the legend
    :type legend_loc: str
    :param legend_fontsize: fontsize of the legends
    :type legend_fontsize: int
    :param marker_size: size of markers
    :type marker_size: int
    :param facecolor: a bool for facecolor of markers
    :type facecolor: bool
    :param subtext: dictionary for subplots text
    :type subtext: dict
    :param canvas: canvas of the figure
    :type canvas: FigureCanvas
    """

    # Create a hidden root window to get screen dimensions
    if use_tk:
        root = tk.Tk()
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        root.destroy()
    else:
        screen_width = 1500
        screen_height = 950


    # Define a scaling factor for the figure size (e.g., 0.8 for 80% of the screen size)
    scaling_factor = 0.8

    # Calculate the figure size based on the screen dimensions and scaling factor
    figure_width = screen_width * scaling_factor
    figure_height = screen_height * scaling_factor

    if new_fig: fig = plt.figure(figsize=(figure_width/100, figure_height/100))
    ax1 = plt.subplot(411)
    plt.plot(s.ppg, 'k', label=None)
    ax2 = plt.subplot(412, sharex=ax1)
    plt.plot(s.vpg, 'k', label=None)
    ax3 = plt.subplot(413, sharex=ax2)
    plt.plot(s.apg, 'k', label=None)
    ax4 = plt.subplot(414, sharex=ax3)
    plt.plot(s.jpg, 'k', label=None)
    if new_fig:
        fig.subplots_adjust(hspace=0, wspace=0)
        canvas = FigureCanvas(fig)

    if len(marker)==0:
        marker = ['o', 's', 's','o', 'o', 's', 'o', 'o', 's', 'o', 's', 'o', 's', 'o', 's']

    color = ['r', 'b', 'g','m', 'r', 'b', 'g', 'r', 'b', 'g', 'm', 'c', 'k', 'r', 'b']

    fid_names = ('sp', 'on', 'dn','dp', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')
    ylabe_names = ['PPG', 'PPG\'', 'PPG\'\'', 'PPG\'\'\'']
    s_type = ['ppg', 'ppg', 'ppg','ppg', 'vpg', 'vpg', 'vpg', 'apg', 'apg', 'apg', 'apg', 'apg', 'apg', 'jpg','jpg']

    str_sig = 0
    end_sig = len(s.ppg)
    len_sig=end_sig-str_sig
    step_small = 1
    if len_sig/s.fs<10:
        step_big = step_small * 1
    else:
        step_big = step_small * 5

    len_sig_sec=int(len_sig / s.fs)
    if len_sig_sec<1:
        len_sig_sec=1

    major_ticks_names = range(0,len_sig_sec,step_big)
    len_ticks_names=len(major_ticks_names)
    major_diff=len_sig/len_ticks_names
    minor_diff = len_sig / len_ticks_names / step_big
    major_ticks = np.arange(str_sig, end_sig, major_diff)
    major_ticks = major_ticks[0:len_ticks_names]
    minor_ticks = np.arange(str_sig, end_sig, minor_diff)

    sig_names=('ppg','vpg','apg','jpg')
    for n in fid_names:
        ind = fid_names.index(n)
        tmp_ind = sig_names.index(s_type[ind])
        if tmp_ind == 0:
            plt_num = 1
            plt.subplot(411)
            plt.title(s.name, fontsize=20)
        else:
            plt_num = tmp_ind+1

        ax = plt.subplot(4, 1, plt_num)

        tmp_pnt=eval("fp." + n + ".values")
        tmp_pnt=tmp_pnt[~np.isnan(tmp_pnt)].astype(int)
        tmp_sig=eval("s." + s_type[ind])

        if facecolor:
            fc=color[ind]
        else:
            fc='none'

        exec("plt.scatter(tmp_pnt, tmp_sig[tmp_pnt], s=marker_size,linewidth=2, marker = marker[ind], facecolors=fc, color=color[ind], label=n)")

        if ind<len(s_type)-1:
            legend_flag = s_type[ind] != s_type[ind + 1]
        else:
            legend_flag=True

        if legend_flag:
            leg = ax.legend(loc=legend_loc, fontsize=legend_fontsize, ncol=2, facecolor="orange", frameon=True, title=title)
            ax.add_artist(leg)

            try:
                tmp_txt=subtext[s_type[ind]]
                plt.text(-.14, 0.95, tmp_txt, fontsize=legend_fontsize,transform=ax.transAxes, va='top', ha='left',
                         bbox={'facecolor': 'yellow', 'alpha': 0.5, 'edgecolor': 'black', 'pad': 5})
            except:
                pass

            if s_type[ind]=='jpg':
                try:
                    tmp_txt = 'MAE(SD):'+'\n'+str(subtext['mae'])+'('+str(subtext['std'])+')'
                    plt.text(-.14, 0, tmp_txt, fontsize=legend_fontsize, transform=ax.transAxes, va='top', ha='left', weight='bold', color='w',
                             bbox={'facecolor': 'red', 'alpha': 0.5, 'edgecolor': 'black', 'pad': 5})
                except:
                    pass

            plt.ylabel(ylabe_names[plt_num - 1], fontsize=20)

            exec("plt.xlim([str_sig,end_sig])")
            for text in leg.get_texts():
                text.set_weight('bold')

            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)

            ax.grid(which='both')
            ax.grid(which='minor', alpha=0.2,axis='x')
            ax.grid(which='major', alpha=1,axis='x')

            plt.yticks([])

    plt.xlabel('Time [s]', fontsize=20)
    plt.xticks(major_ticks, major_ticks_names, fontsize=20)
    if show_fig: plt.show()

    if not(':' in savingfolder):
        relative_path=r'.'+os.sep
    else:
        relative_path=''

    if savefig:
        tmp_dir=savingfolder+os.sep+'PPG_Figures'+os.sep

        os.makedirs(tmp_dir, exist_ok=True)

        canvas.print_png((tmp_dir+'%s_btwn_%s-%s.png') % (s.name,str_sig,end_sig))
        if print_flag: print('Figure has been saved in the "'+tmp_dir+'".')

    return canvas

###########################################################################
################################# Save Data ###############################
###########################################################################

def save_data(savingformat: str, savingfolder: str, print_flag=True, s={}, fp=pd.DataFrame(), bm=pd.DataFrame()):
    """
    Save the results of the filtered PPG analysis.

    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str
    :param savingfolder: location of the saved data
    :type savingfolder: str
    :param print_flag: a bool for print message
    :type print_flag: bool
    :param s: a struct of PPG signal
    :type s: pyPPG.PPG object
    :type s: DotMap
    :param fp: a struct of fiducial points
    :type fp: pyPPG.Fiducial object
    :param bm: a dictionary of biomarkers
    :type bm: pyPPG.Biomarkers object

    :return: file_names: dictionary of the saved file names
    """

    savingfolder = savingfolder.replace('/', '\\')

    if not(':' in savingfolder) and savingfolder[0]!='/':
        relative_path=r'.'+os.sep
    else:
        relative_path=''

    tmp_dir = savingfolder
    os.makedirs(tmp_dir, exist_ok=True)

    temp_dirs = ['Fiducial_points', 'Biomarker_vals', 'Biomarker_stats', 'Biomarker_defs', 'PPG_struct', 'Biomarker_defs_and_stats']
    for i in temp_dirs:
        temp_dir = tmp_dir + os.sep + i + os.sep
        os.makedirs(temp_dir, exist_ok=True)

    keys=s.__dict__.keys()
    keys_list = list(keys)
    sc=DotMap()
    for i in keys_list:
        exec('sc.'+i+' = s.'+i)

    file_names = {}
    file_name = (relative_path + tmp_dir + os.sep + temp_dirs[4] + os.sep + s.name + '_data_btwn_%s-%s.mat')%(s.start_sig,s.end_sig)
    file_names ['data_struct_mat']= file_name
    scipy.io.savemat(file_name, sc)

    try:
        BM_keys = bm.bm_vals.keys()
    except:
        BM_keys={}


    if savingformat=="csv" or savingformat=="both":
        file_name = (relative_path+tmp_dir+os.sep+temp_dirs[0]+os.sep+s.name+'_'+'Fiducials_btwn_%s-%s.csv')%(s.start_sig,s.end_sig)
        file_names ['fiducials_csv']= file_name
        tmp_fp = fp.get_fp()
        tmp_fp.index = tmp_fp.index + 1
        tmp_fp.to_csv(file_name)

        for key in BM_keys:
            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[1]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names [key + '_vals_csv']= file_name
            bm.bm_vals[key].index = bm.bm_vals[key].index + 1
            bm.bm_vals[key].to_csv(file_name,index=True,header=True)

            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[2]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names [key + '_stats_csv']= file_name
            bm.bm_stats[key].to_csv(file_name, index=True, header=True)

            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[3]+os.sep+'%s_btwn_%s-%s.csv')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names [key + '_defs_csv']= file_name
            bm.bm_defs[key].index = bm.bm_defs[key].index + 1
            bm.bm_defs[key].to_csv(file_name, index=True, header=True)

    if savingformat=="mat"  or savingformat=="both":
        file_name = (relative_path+tmp_dir+os.sep+temp_dirs[0]+os.sep+s.name+'_'+'Fiducials_btwn_%s-%s.mat')%(s.start_sig,s.end_sig)
        file_names['fiducials_mat']=file_name
        tmp_fp = fp.get_fp()
        tmp_fp.index = tmp_fp.index + 1
        savemat(file_name, {'PPG_fiducials': tmp_fp.to_records(index=True)})

        for key in BM_keys:

            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[1]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names[key+'_vals_mat']=file_name
            tmp_df_vals=bm.bm_vals[key]
            matlab_struct = tmp_df_vals.to_dict(orient='list')
            savemat(file_name, {'PPG_vals': tmp_df_vals.to_records(index=True)})

            if len(bm.bm_stats)>0:
                file_name = (relative_path+tmp_dir+os.sep+temp_dirs[2]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start_sig,s.end_sig)
                file_names[key + '_stats_mat']=file_name
                tmp_df_stat=bm.bm_stats[key]
                savemat(file_name, {'PPG_stats': tmp_df_stat.to_records(index=True)})

            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[3]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names[key + '_defs_mat']=file_name
            tmp_df_defs=bm.bm_defs[key]
            matlab_struct = tmp_df_defs.to_dict(orient='list')
            savemat(file_name, {'PPG_defs': tmp_df_defs.to_records(index=True)})

            file_name = (relative_path+tmp_dir+os.sep+temp_dirs[5]+os.sep+'%s_btwn_%s-%s.mat')%(s.name+'_'+key,s.start_sig,s.end_sig)
            file_names[key + '_defs_stats_mat']=file_name
            tmp_df_defs2 = tmp_df_defs.drop('name', axis=1)
            tmp_df_defs2.index =tmp_df_defs['name']
            tmp_df_stat2=tmp_df_stat.transpose()
            tmp_df_defs_and_stats=pd.concat([tmp_df_defs2, tmp_df_stat2], axis=1)
            savemat(file_name, {key: tmp_df_defs_and_stats.to_records(index=True)})

    if savingformat != "csv" and savingformat != "mat" and savingformat != "both" and savingformat!="none":
        raise('The file format is not suported for data saving! You can use "mat" or "csv" file formats.')

    if print_flag: print('Results have been saved into the "'+tmp_dir+'".')

    return file_names


def load_fiducials(saved_fiducials=""):
    """
    :param saved_fiducials: path of the matlab struct of the saved fiducial points, where the name field is 'PPG_fiducials'
    :type saved_fiducials: str

    :return: fiducial points: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points
    """
    try:
        loaded_fp = scipy.io.loadmat(saved_fiducials)['PPG_fiducials']
        python_dict = {}
        for field in loaded_fp.dtype.names:
            python_dict[field] = np.squeeze(np.squeeze(loaded_fp[field]))

        fiducials = pd.DataFrame(python_dict)
        for fp_key in python_dict.keys():
            fiducials[fp_key] = [x[0, 0] for x in fiducials[fp_key]]

        return fiducials
    except:
        print('Invalid field name. The supported MATLAB field name for fiducial points is "PPG_fiducials".')
        return

