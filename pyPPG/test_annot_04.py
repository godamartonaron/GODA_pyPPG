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
    ppg_sig_dir='D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/PETE_Matlab'
    ppg_file = '/PPG-BP1.mat'
    annot_path = ppg_sig_dir+'/ANNOTS/MG_PPG-BP_annot/merged'
    sig_path=(ppg_sig_dir + ppg_file)
    input_sig = scipy.io.loadmat(sig_path)

    # Define output variables
    OutData = {}
    set_len=input_sig['ppg_data'].size

    fid_names = ('pk','os','dn', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f')
    dist_error= pd.DataFrame()
    for n in fid_names[2:]:
        exec(n+"=[]")
        exec(n + "r=[]")
        exec(n + "_dist=[]")
        temp_v=np.empty(set_len)
        temp_v[:] = np.NaN
        dist_error[n]=temp_v

    # Flag for plotting: 0 is off, 1 is on
    plt_sig=1

    for i in range(0,set_len):

        # Define sampling frequency, load filtered signal, 1st and 2nd derivative
        fs = input_sig['ppg_data']['fs'][0,0][0][0]
        fs = np.squeeze(fs)

        ppg_v = input_sig['ppg_data']['filt_sig'][0,i]
        ppg_v =np.squeeze(ppg_v)

        drt1 = input_sig['ppg_data']['d1'][0,i]
        drt1 = np.squeeze(drt1)

        drt2 = input_sig['ppg_data']['d2'][0,i]
        drt2 = np.squeeze(drt2)

        # Plot filtered signal, 1st and 2nd derivative
        if plt_sig==1:
            fig = plt.figure(figsize=(15, 7))
            plt.plot(ppg_v,'r',label='x')
            plt.plot(drt1,'b',label='dx')
            plt.plot(drt2,'k',label='ddx')

        # Load annotated fiducial points
        name=input_sig['ppg_data']['name'][0, i][0]
        annot_file=annot_path+'/'+name+'.mat'
        annot=scipy.io.loadmat(annot_file)

        for n in fid_names:
            exec("ref_" + n +" = np.squeeze(np.round(annot['annot']['"+n+"'][0, 0][0][0]['t'] * fs).astype(int))")

        # Plot onset and peak
        pks,ons=[],[]
        exec("pks = [ref_pk]")
        exec("ons = ref_os")
        if plt_sig==1:
            plt.scatter(pks, ppg_v[pks], s=60, linewidth=2, marker='o',  facecolors='c', edgecolors='r', label='pk')
            plt.scatter(ons, ppg_v[ons], s=60, linewidth=2, marker='s',  facecolors='c', edgecolors='b', label='os')

        # Detect fiducial points
        det_dn=getDicroticNotch(ppg_v, fs, pks, ons)
        det_u, det_v, det_w = getFirstDerivitivePoints(ppg_v, fs, ons)
        det_a, det_b, det_c, det_d, det_e, det_f = getSecondDerivitivePoints(ppg_v,fs, ons)

        for n in fid_names[2:]:
            # Calculate distance error
            exec (n+".append(np.squeeze(det_"+n+"))")
            exec("temp_dist=np.squeeze(ref_"+n+" - det_"+n+")")
            if eval("temp_dist.size") > 0:
                exec(n + "_dist.append(temp_dist)")
                exec ("dist_error['"+n+"'][i] = temp_dist")

            # Plot fiducial points
            ind = fid_names.index(n)-2
            s_type = ['ppg_v', 'drt1', 'drt1', 'drt1', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2']
            marker = ['s', 'x','o', '*']*5
            color = ['b', 'r', 'c', 'm', 'k', 'g', 'm', 'b', 'r', 'c', 'm']

            is_fidu=0
            exec("is_fidu=~np.isnan(np.squeeze(det_" + n + "))")

            if  plt_sig==1 and is_fidu:
                exec("plt.scatter(ref_" + n + "," + s_type[ind] + "[ref_" + n + "], s=60,linewidth=2, marker = marker[ind*2], facecolors='none', edgecolors=color[ind], label='ref " + n + "')")
                exec("plt.scatter(det_" + n + "," + s_type[ind] + "[det_" + n + "], s=60,linewidth=2, marker = marker[ind*2+1], color=color[ind+1], label='det " + n + "')")

        # Plot show
        if plt_sig == 1:
            plt.legend(loc=4, prop={'size': 10})
            plt.title(name, fontsize=20)
            plt.xlabel('Time [ms]', fontsize=20)
            plt.ylabel('Pulse Wave', fontsize=20)
            plt.grid(color='g', linestyle='--', linewidth=0.5)
            plt.savefig(('temp_dir/figs/py_%s.png')%(name))
            # plt.show()
            plt.close('all')

        # Print distance error
        print('\n',name,' DIFF:')
        test_list=dist_error.iloc[i].values
        if sum(np.isnan(test_list))==0:
            temp_error=dist_error.iloc[i].astype(int)
            df = pd.DataFrame(temp_error).transpose()
            print(df)

        # Save distance error in .mat file
        file_name = 'temp_dir/PPG-BP1_eval02.mat'
        OutData['dist_error'] = dist_error.to_numpy()
        scipy.io.savemat(file_name, OutData)

    # Calculate Mean Absolute Error (MAE), Standard Deviation (STD), and Mean Error (BIAS)
    MAE={}
    STD={}
    BIAS={}
    for n in fid_names[2:]:
        exec("MAE['" + n + "'] = np.round(np.nanmean(np.absolute(" + n + "_dist)),2)")
        exec("STD['" + n + "'] = np.round(np.nanstd(" + n + "_dist), 2)")
        exec("BIAS['" + n + "'] = np.round(np.nanmean(" + n + "_dist),2)")

    # Print results
    print('-------------------------------------------')
    print('MAE',': ', MAE)
    print('STD',': ', STD)
    print('BIAS', ': ', BIAS)


