from pyPPG.example import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap

from scipy.signal import filtfilt, find_peaks

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
    cmp_pc = 0

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
        print(i)

        # Define sampling frequency, load filtered signal, 1st and 2nd derivative
        fs = input_sig['ppg_data']['fs'][0,0][0][0]
        fs = np.squeeze(fs)

        ppg_v = input_sig['ppg_data']['filt_sig'][0,i]
        ppg_v = input_sig['ppg_data']['sig'][0, i]
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
        if plt_sig==1:
            plt.scatter(pks, ppg_v[pks], s=150, linewidth=2, marker='o',  facecolors='c', edgecolors='r', label='pk')
            plt.scatter(ons, ppg_v[ons], s=150, linewidth=2, marker='s',  facecolors='m', edgecolors='b', label='on')

        # Detect fiducial points
        s = DotMap()
        s.fs = fs
        s.v = ppg_v
        s.ppg = ppg_v
        s.vpg = drt1
        s.apg = drt2
        s.jpg = drt3
        s.name = name

        # s.filt_d1 = drt1
        # s.filt_d2 = drt2
        # s.filt_d3 = drt3

        if cmp_pc !=1:

            ## Preprocessing
            sm_wins = {'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10}
            #prep = PP.Preprocess(fL=0, fH=12, order=4, sm_wins=sm_wins)
            prep = PP.Preprocess(fL=0.5, fH=12, order=4, sm_wins=sm_wins)
            prep = PP.Preprocess(fL=0.5000001, fH=12, order=4, sm_wins=sm_wins)
            #prep = PP.Preprocess(fL=0.49999999, fH=12, order=4, sm_wins=sm_wins)

            # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
            s.filtering = True
            # s.v = np.concatenate((ppg_v,ppg_v+(ppg_v[-1]-ppg_v[0])),axis=0)
            # s.name = name+'+'+name
            s.ppg, s.vpg, s.apg, s.jpg = prep.get_signals(s)

            # drt0 = input_sig['ppg_data']['filt_sig'][0, i]
            # plt.plot(s.ppg / np.max(s.ppg), label='p.ppg')
            # plt.plot(drt0 / np.max(drt0), label='m.ppg')
            # plt.legend()
            # plt.show()

            ## Create a PPG class

            s.correct = True#False
            s_class = PPG(s,check_ppg=False)


            FP.FpCollection(s_class)
            fpex = FP.FpCollection(s_class)

            # Extract fiducial points
            det_dn = np.array(fpex.get_dicrotic_notch(pks, ons))
            drt1_fp = fpex.get_vpg_fiducials(ons)
            drt2_fp = fpex.get_apg_fiducials(ons, pks)
            drt3_fp = fpex.get_jpg_fiducials(ons,drt2_fp)

            det_dp = fpex.get_diastolic_peak(ons, det_dn, drt2_fp.e)

            ## Added 10/09/2023
            drt1_fp[drt1_fp.columns[drt1_fp.isna().any()]] = 0
            drt2_fp[drt2_fp.columns[drt2_fp.isna().any()]] = 0
            drt3_fp[drt3_fp.columns[drt3_fp.isna().any()]] = 0
            type_sig=['drt1_fp']*3+['drt2_fp']*6+['drt3_fp']*2

            for j in range(3,len(fid_names)):
                try:
                    exec("det_"+fid_names[j]+"=np.array(int("+type_sig[j-3]+"."+fid_names[j]+"))")
                except:
                    exec("det_" + fid_names[j] + "=np.NaN")

            ### Check fidu
            if det_a>75:
                win_on=75
            else:
                win_on = det_a

            det_on = np.argmax(drt3[det_a-win_on:det_a])+det_a-win_on

            ## Un comment 10/09/2023
            strt_dn=det_e
            stp_dn=det_f

            if det_w>det_f:
                det_w = det_f

            if det_w<det_e:
                det_w = np.argmax(drt1[det_e:det_f])+det_e
            ##


            ## Comment 10/09/2023
            # # Correct v-point
            if det_v>det_e:
                try:
                    det_v = np.argmin(drt1[det_u:det_e])+ det_u
                    det_w = find_peaks(drt1[det_v:det_f])[0][0] + det_v
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

        elif cmp_pc==1:
            for n in f_names:
                if n!='on':
                    exec("M_FID_2['" + n + "'][i] =det_" + n)
                else:
                    exec("M_FID_2['" + n + "'][i] =np.mean(det_" + n+")")

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

                # fid_ref= pd.DataFrame()
                # for fi in fid_names:
                #     try:
                #         fid_ref[fi] = eval("ref_"+fi)
                #     except:
                #         fid_ref[fi] = np.nan
                #
                # fid_ref['dp']=np.nan
                # fp = Fiducials(fp=fid_ref)
                # marker1 = ['o', 's', 's','o', 'o', 's', 'o', 'o', 's', 'o', 's', 'o', 's', 'o', 's']
                # plot_fiducials(s=s, fp=fp, savefig=False, show_fig=False, print_flag=False, use_tk=False)
                #
                # fid_det = pd.DataFrame()
                # for fi in fid_names:
                #     try:
                #         fid_det[fi] = eval("det_"+fi)
                #     except:
                #         fid_det[fi] = np.nan
                #
                # fid_det['dp']=np.nan
                # fp = Fiducials(fp=fid_det)
                # marker2 = ['*', 'x', 'x', '*','*', 'x', '*', '*', 'x', '*', 'x', '*', 'x', '*', 'x']
                # plot_fiducials(s=s, fp=fp, savefig=False, show_fig=False, print_flag=False, use_tk=False, new_fig=False, marker=marker2)
                #
                # plt.show()

                # s_type = ['s.ppg']*2+['s.vpg']*3+['s.apg']*6+['j.apg']*2
                # plt_nums = [1]*2+[2]*3+[3]*6+[4]*2

                if  plt_sig==1 and is_fidu:
                    if s_type[ind][-1]=='v':
                        plt_num=1
                        plt.subplot(411)
                        plt.title(name, fontsize=20)
                    else:
                        plt_num=int(s_type[ind][-1])+1
                        # plt_num = plt_nums[ind]

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
                file_name = 'temp_dir/results/MG-PC_errors.mat'
            else:
                OutData['pyPPG_FID'] = M_FID_2.to_numpy()
                try:
                    if len(annot_path2)>0:
                        file_name = 'temp_dir/results/MG-pyPPG_errors2.mat'
                except:
                    file_name = 'temp_dir/results/PC-pyPPG_errors2.mat'
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
