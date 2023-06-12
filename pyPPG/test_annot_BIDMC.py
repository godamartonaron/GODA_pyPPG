from FiducialPoints import*
from Biomarkers import*
from Statistics import*
from Summary import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time

from six.moves import cPickle as pickle

import matplotlib.mlab
from scipy.io import savemat

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':

    # Define input directories and files
    ppg_sig_dir = 'D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/GIT_PPG_annot'
    ppg_file = '/BIDMC-1min.mat'
    annot_path = ppg_sig_dir + '/temp_dir/BIDMC'

    sig_path=(ppg_sig_dir + ppg_file)
    input_sig = scipy.io.loadmat(sig_path)

    # Define output variables
    OutData = {}

    # Flag for plotting: 0 is off, 1 is on
    plt_sig = 1

    # Flag for comparison with PC: 0 is off, 1 is on
    cmp_pc = 1

    fid_names = ('pk', 'os')#, 'dn', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e')#, 'f','p1','p2')
    if cmp_pc!=1:
        f_names=fid_names [2:]
    else:
        f_names = fid_names[:]

    STAT_MX = {}
    for n in f_names:
        exec(n+"=[]")
        exec(n + "r=[]")
        exec("dist_"+ n +"=[]")

        STAT_MX[n] = {}
        STAT_MX[n]["MAE"] = {}
        STAT_MX[n]["STD"] = {}
        STAT_MX[n]["MBE"] = {}

    rec_set=np.array([18,23,25])-1
    for i in rec_set:

        # Define sampling frequency, load filtered signal, 1st and 2nd derivative
        fs = input_sig['ppg_data']['fs'][0,0][0][0]
        fs = np.squeeze(fs)

        ppg_v = input_sig['ppg_data']['filt_sig'][0,i]
        ppg_v =np.squeeze(ppg_v)

        drt1 = input_sig['ppg_data']['d1'][0,i]
        drt1 = np.squeeze(drt1)

        drt2 = input_sig['ppg_data']['d2'][0,i]
        drt2 = np.squeeze(drt2)

        win = (fs * 0.01).astype(int)
        B = 1 / win * np.ones(win)
        drt3 = np.gradient(drt2)
        # drt3 = filtfilt(B, 1, drt3)
        drt3 = drt3 / max(drt3) * max(ppg_v)

        # Plot filtered signal, 1st and 2nd derivative
        if plt_sig==1:
            fig = plt.figure(figsize=(15, 7))
            plt.plot(ppg_v,'k',label='PPG')
            # plt.plot(drt1,'k--',label='PPG''')
            # plt.plot(drt2,'k:',label='PPG"')
            # plt.plot(drt3, 'k-.', label='PPG''"')

        # Load annotated fiducial points
        name=input_sig['ppg_data']['name'][0, i][0]
        annot_file=annot_path+'/'+name+'.mat'
        annot=scipy.io.loadmat(annot_file)

        for n in fid_names:
            exec("ref_" + n +" = np.squeeze(np.round(annot['annot']['"+n+"'][0, 0][0][0]['t'] * fs).astype(int))")

        if cmp_pc==1:
            s = DotMap()
            s.fs = fs
            s.filt_sig = ppg_v
            s.filt_d1 = drt1
            s.filt_d2 = drt2
            s.filt_d3 = drt3
            peak_detector = 'abp'
            det_pk, det_os = abdp_beat_detector(s, peak_detector)
            det_os = np.array(det_os)

        # Plot onset and peak
        pks,ons=[],[]
        exec("pks = [ref_pk]")
        exec("ons = ref_os")
        if plt_sig==1:
            plt.scatter(pks, ppg_v[pks], s=150, linewidth=2, marker='o', facecolors='w', edgecolors='r', label='ref pk')
            plt.scatter(ons, ppg_v[ons], s=150, linewidth=2, marker='s', facecolors='w', edgecolors='b', label='ref os')
            plt.scatter(det_pk, ppg_v[det_pk], s=150, linewidth=2, marker='x', color='g', label='det pk')
            plt.scatter(det_os, ppg_v[det_os], s=150, linewidth=2, marker='+', color='m', label='det os')

        print('\n', name, ' DIFF:')
        for n in f_names:
            # Calculate distance error
            if eval("det_"+n+".size") > 0:
                exec("len_n=len(det_" + n + ")")
                for j in range(len_n):
                    exec("temp_dist=min(abs(det_"+n+"[j]-ref_"+n+"))")
                    exec("dist_" + n + ".append(temp_dist)")

            exec("MAE= np.round(np.nanmean(np.absolute(dist_" + n + ")),2)")
            exec("STD = np.round(np.nanstd(dist_" + n + "), 2)")
            exec("MBE= np.round(np.nanmean(dist_" + n + "),2)")

            STAT_MX[n]["MAE"][name] = MAE
            STAT_MX[n]["STD"][name] = STD
            STAT_MX[n]["MBE"][name] = MBE

        # Plot show
        if plt_sig == 1:
            plt.legend(loc=4, prop={'size': 11.8})
            plt.title(name, fontsize=20)
            plt.xlabel('Time [ms]', fontsize=20)
            plt.ylabel('Amplitude [nu]', fontsize=20)
            plt.yticks([])
            plt.xticks(fontsize=20)
            plt.show()
            plt.close('all')

    for n in f_names:
        STAT_MX[n] = pd.DataFrame.from_dict(STAT_MX[n])
        print(n)
        print(STAT_MX[n])

    print('Program finished!')