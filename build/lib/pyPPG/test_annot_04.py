from FiducialPoints import*
from Biomarkers import*
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
from scipy.io import savemat

###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':

    # Define input directories and files
    ppg_sig_dir = 'D:/ALL_DATA/Uni/Subjects/ITK_Adjunktus/HAIFA/TECHNION-BME/Research/PPG/GIT_PPG_annot'
    ppg_file = '/PPG-BP1.mat'
    annot_path = ppg_sig_dir + '/temp_dir/MG05_PPG-BP_annot/merged' #'/ANNOTS/MG_PPG-BP_annot/merged'
    annot_path2 = ppg_sig_dir + '/temp_dir/PC05_PPG-BP_annot/merged'
    sig_path=(ppg_sig_dir + ppg_file)
    input_sig = scipy.io.loadmat(sig_path)

    # Define output variables
    OutData = {}
    set_len=input_sig['ppg_data'].size

    # Flag for plotting: 0 is off, 1 is on
    plt_sig = 0

    # Flag for comparison with PC: 0 is off, 1 is on
    cmp_pc = 1

    # Flag for saving: 0 is no, 1 is yes
    save = 1

    fid_names = ('sp', 'on', 'dn', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f','p1','p2')
    if cmp_pc:
        f_names=fid_names [:]
    else:
        f_names = fid_names[1:]

    dist_error= pd.DataFrame()
    M_FID_1 = pd.DataFrame()
    M_FID_2 = pd.DataFrame()
    for n in f_names:
        exec(n+"=[]")
        exec(n + "r=[]")
        exec("dist_"+ n +"=[]")
        temp_v=np.empty(set_len)
        temp_v[:] = np.NaN
        dist_error[n]=temp_v
        M_FID_1[n] = temp_v
        M_FID_2[n] = temp_v

    annot_error = pd.DataFrame()
    for n in fid_names:
        annot_error[n] = temp_v

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

        win = (fs * 0.01).astype(int)
        B = 1 / win * np.ones(win)
        drt3 = np.gradient(drt2)
        drt3 = filtfilt(B, 1, drt3)
        drt3 = drt3 / max(drt3) * max(ppg_v)

        # Plot filtered signal, 1st and 2nd derivative
        if plt_sig:
            fig = plt.figure(figsize=(8, 9))
            ax1 = plt.subplot(411)
            plt.plot(ppg_v,'k',label=None)
            ax2 = plt.subplot(412,sharex=ax1)
            plt.plot(drt1,'k--',label=None)
            ax3 = plt.subplot(413,sharex=ax2)
            plt.plot(drt2,'k:',label=None)
            ax4 = plt.subplot(414,sharex=ax3)
            plt.plot(drt3, 'k-.', label=None)
            fig.subplots_adjust(hspace=0, wspace=0)

        # Load annotated fiducial points
        name=input_sig['ppg_data']['name'][0, i][0]
        annot_file=annot_path+'/'+name+'.mat'
        annot=scipy.io.loadmat(annot_file)

        for n in fid_names:
            if n=='sp':
                exec("ref_" + n + " = np.array(np.squeeze(np.round(annot['annot']['pk'][0, 0][0][0]['t'] * fs).astype(int)))")
            elif n=='on':
                exec("ref_" + n + " = np.array(np.squeeze(np.round(annot['annot']['os'][0, 0][0][0]['t'] * fs).astype(int)))")
            else:
                exec("ref_" + n +" = np.array(np.squeeze(np.round(annot['annot']['"+n+"'][0, 0][0][0]['t'] * fs).astype(int)))")

            if eval("ref_" + n + ".size") > 1 and n != 'on':
                exec("ref_" + n + "=  np.array(ref_" + n + "[0])")
                # exec("ref_" + n + "=  np.array([ref_" + n + "[0]])")
                exec("annot_error['" + n + "'][i] =1")
                print(n + ' error: ', name)
            elif eval("ref_" + n + ".size") < 1 and n != 'on':
                exec("ref_" + n + "= np.array(np.NaN)")
                exec("annot_error['" + n + "'][i] =1")
                print(n + ' error: ', name)

        if ref_on.size > 2:
            exec("annot_error['on'][i] =1")
            print('on error: ', name)
            exec("ref_on = ref_on[0:2]")

        if ref_on.size < 2:
            exec("annot_error['on'][i] =1")
            print('on error: ', name)
            exec("ref_on = np.squeeze([ref_on,ref_on])")

        for n in f_names:
            if n=="on":
                exec("M_FID_1['" + n + "'][i] =ref_" + n+"[0]")
            else:
                exec("M_FID_1['" + n + "'][i] =ref_" + n)

        if cmp_pc:
            annot_file2=annot_path2+'/'+name+'.mat'
            annot2=scipy.io.loadmat(annot_file2)

            for n in fid_names:
                if n == 'sp':
                    exec("det_" + n + " = np.array(np.squeeze(np.round(annot['annot']['pk'][0, 0][0][0]['t'] * fs).astype(int)))")
                elif n == 'on':
                    exec("det_" + n + " = np.array(np.squeeze(np.round(annot['annot']['os'][0, 0][0][0]['t'] * fs).astype(int)))")
                else:
                    exec("det_" + n +" = np.squeeze(np.round(annot2['annot']['"+n+"'][0, 0][0][0]['t'] * fs).astype(int))")

                if eval("det_" + n +".size")>1 and n!='on':
                    exec("det_"+n+"= np.array(det_"+n+"[0])")
                    exec ("annot_error['"+n+"'][i] =1")
                    print(n+' error: ', name)
                elif eval("det_" + n +".size")<1 and n!='on':
                    exec("det_"+n+"= np.NaN")
                    exec ("annot_error['"+n+"'][i] =1")
                    print(n+' error: ', name)

            if det_on.size > 2:
                exec("annot_error['on'][i] =1")
                print('on error: ', name)
                exec("det_on = det_on[0:2]")

            if det_on.size < 2:
                exec("annot_error['on'][i] =1")
                print('on error: ', name)
                exec("det_on = np.squeeze([det_on,det_on])")

            for n in f_names:
                if n == "on":
                    exec("M_FID_2['" + n + "'][i] =det_" + n + "[0]")
                else:
                    exec("M_FID_2['" + n + "'][i] =det_" + n)

        # Plot onset and peak
        pks,ons=[],[]
        exec("pks = [ref_sp]")
        exec("ons = ref_on")
        # if plt_sig==1:
        #     plt.scatter(pks, ppg_v[pks], s=150, linewidth=2, marker='o',  facecolors='c', edgecolors='r', label='pk')
        #     plt.scatter(ons, ppg_v[ons], s=150, linewidth=2, marker='s',  facecolors='m', edgecolors='b', label='on')

        # Detect fiducial points
        s = DotMap()
        s.fs = fs
        s.filt_sig = ppg_v
        s.filt_d1 = drt1
        s.filt_d2 = drt2
        s.filt_d3 = drt3

        if cmp_pc !=1:

            det_dn = np.array(getDicroticNotch(s, pks, ons))
            drt1_fp = getFirstDerivitivePoints(s, ons)
            drt2_fp = getSecondDerivitivePoints(s, ons, pks)
            drt3_fp = getThirdDerivitivePoints(s, ons,drt2_fp)

            det_dp = getDiastolicPeak(s, ons, det_dn, drt2_fp.e)

            det_u, det_v, det_w = np.array(int(drt1_fp.u)), np.array(int(drt1_fp.v)), np.array(int(drt1_fp.w))
            det_a, det_b, det_c, det_d, det_e, det_f = np.array(int(drt2_fp.a)), np.array(int(drt2_fp.b)), np.array(int(drt2_fp.c)), np.array(int(drt2_fp.d)), np.array(int(drt2_fp.e)), np.array(int(drt2_fp.f))
            det_p1, det_p2 = np.array(int(drt3_fp.p1)), np.array(int(drt3_fp.p2))

            ### Check fidu
            if det_a>75:
                win_on=75
            else:
                win_on = det_a

            det_on = np.argmax(drt3[det_a-win_on:det_a])+det_a-win_on

            # strt_dn=det_e
            # stp_dn=det_f
            # det_dn = np.argmin(drt3[strt_dn:stp_dn])+strt_dn

            # if det_w>det_f:
            #     det_w = det_f
            #
            # if det_w<det_e:
            #     det_w = np.argmax(drt1[det_e:det_f])+det_e


            try:
                temp_segment = s.filt_sig[int(ref_sp):int(det_dp)]
                min_dn = find_peaks(-temp_segment)[0] + ref_sp
                diff_dn = abs(min_dn - det_dp)
                if len(min_dn) > 0 and diff_dn > round(s.fs / 100):
                    try:
                        strt_dn = int(ref_sp)
                        stp_dn = int(det_f)
                        det_dn = find_peaks(-s.filt_sig[strt_dn:stp_dn])[0][-1] + strt_dn
                        if det_dn > min_dn:
                            det_dn = min_dn
                    except:
                        strt_dn = det_e
                        stp_dn = det_f
                        det_dn = np.argmin(drt3[strt_dn:stp_dn])+strt_dn
                        if det_dn > min_dn:
                            det_dn = min_dn
            except:
                pass

            # Correct v-point
            if det_v>det_e:
                det_v = np.argmin(drt1[det_u:det_e])+ det_u
                det_w = find_peaks(drt1[det_v:det_f])[0][0] + det_v

            # Correct w-point
            try:
                temp_end = int(np.diff(ons) * 0.8)
                temp_segment = s.filt_d1[int(det_dn):int(ons[0] + temp_end)]
                min_w = find_peaks(-temp_segment)[0] + det_dn
                if len(min_w)>1:
                    min_w=min_w[0]

                if det_w<det_e:
                    det_w = np.argmax(drt1[det_e:det_f])+det_e

                if det_w>det_f:
                    det_w = det_f

                if det_w > min_w:
                    det_w = min_w
            except:
                pass

            # Correct f-point
            try:
                temp_end = int(np.diff(ons) * 0.8)
                temp_segment = s.filt_d2[int(det_e):int(ons[0] + temp_end)]
                min_f = np.argmin(temp_segment) + det_e

                if det_w > det_f:
                    det_f = min_f
            except:
                pass

            for n in f_names:
                exec("M_FID_2['" + n + "'][i] =det_" + n)

            # Calculate distance error
            for n in f_names:
                exec (n+".append(np.squeeze(det_"+n+"))")
                if eval("ref_"+n+".size") > 0:
                    if n == 'on':
                        if  cmp_pc:
                            exec("temp_dist=np.array(np.mean(np.squeeze(ref_"+n+" - det_"+n+")))")
                        else:
                            exec("temp_dist=np.array(np.mean(np.squeeze(ref_" + n + "[0] - det_" + n + ")))")
                    else:
                        exec("temp_dist=np.squeeze(ref_"+n+" - det_"+n+")")
                    if eval("temp_dist.size") > 0:
                        exec("dist_"+ n + ".append(temp_dist)")
                        exec ("dist_error['"+n+"'][i] = temp_dist")

        # Plot fiducial points
        if cmp_pc:
            str_i = 0
        else:
            str_i = 1

        plt_typ = ["ref", "det"]

        zero_sig=fs/10
        if eval("zero_sig<ref_on[0]"):
            str_sig=eval("int(ref_on[0]-zero_sig)")
        else:
            str_sig = eval("int(ref_on[0])")

        end_sig=eval("int(ref_on[0] + 1.5 * np.diff(ref_on))")

        str_sig=0
        end_sig=len(ppg_v)
        len_sig=end_sig-str_sig

        # ref_on = ref_on[0]

        for m in plt_typ:
            for n in f_names:
                ind=f_names.index(n)
                if cmp_pc != 1:
                    s_type = ['ppg_v', 'ppg_v', 'drt1', 'drt1', 'drt1', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2', 'drt3', 'drt3']
                else:
                    s_type = ['ppg_v', 'ppg_v', 'ppg_v', 'drt1', 'drt1', 'drt1', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2', 'drt2','drt3', 'drt3']

                label = m

                # if m=='det':
                #     label='ref2'
                # else:
                #     label = 'ref1'

                # marker = ['s', 'x', 'o', '+']*7
                # color = ['b', 'r', 'c', 'm', 'k', 'g', 'm', 'b', 'r', 'c', 'g', 'k', 'b','c','r']

                marker1 = ['o', 's', 's', 'o', 's', 'o', 'o', 's', 'o', 's', 'o', 's', 'o', 's']
                marker2 = ['*', 'x', 'x', '*', 'x', '*', '*', 'x', '*', 'x', '*', 'x', '*', 'x']
                color = ['r', 'b', 'g', 'r', 'b', 'g', 'r', 'b', 'g', 'm', 'c', 'k', 'r','b']

                if m==plt_typ[0]:
                    marker = marker1[str_i:]
                    color = color[str_i:]
                    facecolors='none'
                    sm = 120
                elif m==plt_typ[1]:
                    marker = marker2[str_i:]
                    color = color[str_i:]
                    facecolors = color[ind]
                    sm = 50

                ylabe_names=['PPG','PPG\'','PPG\'\'','PPG\'\'\'']

                is_fidu=0
                if n == 'on':
                    if cmp_pc:
                        exec("is_fidu=~min(np.isnan(np.squeeze(det_" + n + "))) and ~min(np.isnan(np.squeeze(ref_" + n + ")))")
                    else:
                        exec("is_fidu=~np.isnan(np.squeeze(det_" + n + ")) and ~np.isnan(np.squeeze(ref_" + n + "[0]))")
                else:
                    exec("is_fidu=~np.isnan(np.squeeze(det_" + n + ")) and ~np.isnan(np.squeeze(ref_" + n + "))")

                # exec("is_fidu=~np.isnan(np.squeeze(det_" + n + ")) and ~np.isnan(np.squeeze(ref_" + n + "))")

                if  plt_sig==1 and is_fidu:
                    if s_type[ind][-1]=='v':
                        plt_num=1
                        plt.subplot(411)
                        plt.title(name, fontsize=20)
                    else:
                        plt_num=int(s_type[ind][-1])+1

                    plt.subplot(4,1,plt_num)

                    exec("plt.scatter("+m+"_" + n + "," + s_type[ind] + "["+m+"_" + n + "], s=sm,linewidth=2, marker = marker[ind],facecolors=facecolors, color=color[ind], label= label+' " + n + "')")
                    #exec("plt.scatter("+m+"_" + n + "," + s_type[ind] + "["+m+"_" + n + "], s=sm,linewidth=2, marker = marker[ind],facecolors=facecolors, color=color[ind], label=None)")
                    plt.ylabel(ylabe_names[plt_num-1], fontsize=20)
                    plt.yticks([])
                    exec("plt.ylim([min("+s_type[ind]+")-200,max("+s_type[ind]+")+200])")
                    exec("plt.xlim([str_sig,end_sig])")

                    plt.legend(loc='upper right', fontsize=12, ncol=2)

        # Plot show
        if plt_sig == 1:
            plt.xlabel('Time [ms]', fontsize=20)
            plt.xticks(fontsize=20)
            plt.xticks(range(str_sig, end_sig, 500),range(0, len_sig, 500))
            plt.xticks(range(str_sig, end_sig, 100))
            # plt.grid(color='g', linestyle='--', linewidth=0.5)
            # plt.savefig(('temp_dir/MG_PC_annot5/py_%s.png')%(name))
            plt.show()
            plt.close('all')

        # Print distance error
        print('\n',name,' DIFF:')
        test_list=dist_error.iloc[i].values
        if sum(np.isnan(test_list))==0:
            temp_error=dist_error.iloc[i].astype(int)
            df = pd.DataFrame(temp_error).transpose()
            print(df)

        # Save distance error and annotation error in .mat file
        OutData['dist_error'] = dist_error.to_numpy()
        OutData['annot_error'] = annot_error.to_numpy()
        OutData['MG_FID'] = M_FID_1.to_numpy()
        if save:
            if cmp_pc:
                OutData['PC_FID'] = M_FID_2.to_numpy()
                file_name = 'temp_dir/MG-PC_errors.mat'
            else:
                OutData['pyPPG_FID'] = M_FID_2.to_numpy()
                try:
                    if len(annot_path2)>0:
                        file_name = 'temp_dir/MG-pyPPG_errors2.mat'
                except:
                    file_name = 'temp_dir/PC-pyPPG_errors2.mat'
            scipy.io.savemat(file_name, OutData)

    # Calculate Mean Absolute Error (MAE), Standard Deviation (STD), and Mean Error (BIAS)
    MAE={}
    STD={}
    BIAS={}
    for n in f_names:
        exec("MAE['" + n + "'] = np.round(np.nanmean(np.absolute(dist_"+ n + ")),2)")
        exec("STD['" + n + "'] = np.round(np.nanstd(dist_"+ n + "), 2)")
        exec("BIAS['" + n + "'] = np.round(np.nanmean(dist_"+ n + "),2)")

    # Print results
    print('-------------------------------------------')
    print('MAE',': ', MAE)
    print('STD',': ', STD)
    print('BIAS', ': ', BIAS)
    print('Program finished!')

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


