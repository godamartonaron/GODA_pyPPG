import copy
import os

import pandas as pd

import pyPPG
from pyPPG.example import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap

from scipy.signal import filtfilt, find_peaks
from pyPPG.ppg_bm.statistics import get_statistics
from datetime import datetime


class PulseWaveAnal:

    ###########################################################################
    ###################### Initialization Pulse Wave Analysis #################
    ###########################################################################
    def __init__(self):
        """
        The purpose of the PulseWaveAnal class is to:
            - validate the fiducial point detection
            - compare the fiducial point annotation
            - extract the pulse wave biomarkers    

        """


    ###########################################################################
    ###################### Get Reference Fiducial Points ######################
    ###########################################################################
    def get_ref_fp(self, name, fs, fid_names, annot):
        error = {key: 0 for key in fid_names}

        for n in fid_names:
            try:
                if n == 'sp':
                    exec("ref_" + n + " = np.array(np.squeeze(np.round(annot['annot']['pk'][0, 0][0][0]['t'] * fs).astype(int)))")
                elif n == 'on':
                    ref_on = eval("np.array(np.squeeze(np.round(annot['annot']['os'][0, 0][0][0]['t'] * fs).astype(int)))")
                else:
                    exec("ref_" + n + " = np.array(np.squeeze(np.round(annot['annot']['" + n + "'][0, 0][0][0]['t'] * fs).astype(int)))")

                if eval("ref_" + n + ".size") > 1 and n != 'on':
                    exec("ref_" + n + "=  np.array(ref_" + n + "[0])")
                    exec("error['" + n + "'] =1")
                    print('\n'+n + ' error: ', name)
                elif eval("ref_" + n + ".size") < 1 and n != 'on':
                    exec("ref_" + n + "= np.array(np.NaN)")
                    exec("error['" + n + "'] =1")
                    print(n + ' error: ', name)
            except:
                pass

        if ref_on.size > 2:
            exec("error['on']=1")
            print('on error: ', name)
            exec("ref_on = ref_on[0:2]")

        if ref_on.size < 2:
            exec("error['on']=1")
            print('on error: ', name)
            exec("ref_on = np.squeeze([ref_on,ref_on])")

        ref_fp= {}
        for n in fid_names:
            if n=="on":
                exec("ref_fp['on'] = ref_" + n+"[0]")
                exec("ref_fp['off'] = ref_" + n + "[1]")
            else:
                try:
                    exec("ref_fp['" + n + "'] = int(ref_" + n+")")
                except:
                    exec("ref_fp['" + n + "'] = np.nan")

        return ref_fp, error

    ###########################################################################
    ########################## Merge Fiducial Points ##########################
    ###########################################################################
    def merge_fiducials(self, ppg_fp, vpg_fp, apg_fp, jpg_fp):
        fiducials = pd.DataFrame()
        for temp_sig in (ppg_fp, vpg_fp, apg_fp, jpg_fp):
            for key in list(temp_sig.keys()):
                fiducials[key] = temp_sig[key].values

        return fiducials

    ###########################################################################
    ######################### Calculate Distance Error ########################
    ###########################################################################
    def get_dist_error(self, fp1,fp2,compare):
        d_error={key: 0 for key in fp1.columns}

        for n in fp1.columns:
            if eval("fp1." + n + ".size") > 0:
                exec("temp_dist=(fp1." + n + " - fp2." + n + ")[0]")
                if eval("temp_dist.size") > 0:
                    exec("d_error['" + n + "'] = temp_dist")

        return d_error

    ###########################################################################
    ########################### Plot Fiducial Points ##########################
    ###########################################################################

    def plot_pulse_wave(self, s=pd.DataFrame(),fp1=pd.DataFrame(),fp2=pd.DataFrame(),d_error=pd.DataFrame(),compare=False,show_fig=False,dname='',annot1='', annot2=''):
        marker1 = ['o', 's', 's', 'o', 'o', 's', 'o', 'o', 's', 'o', 's', 'o', 's', 'o', 's']
        marker2 = ['*', 'X', 'X', '*', '*', 'X', '*', '*', 'X', '*', 'X', '*', 'X', '*', 'X']
        legend_fontsize = 9

        subtext = {}

        fp_names = {'ppg': ['on', 'sp', 'dn'],
                    'vpg': ['u', 'v', 'w'],
                    'apg': ['a', 'b', 'c', 'd', 'e', 'f'],
                    'jpg': ['p1', 'p2']}

        for tmp_sig in fp_names:
            tmp_txt = 'Dist. errors'
            for tmp_fp in fp_names[tmp_sig]:
                try:
                    tmp_txt = tmp_txt + '\n  ' + tmp_fp + ': ' + str(int(d_error[tmp_fp]))
                except:
                    tmp_txt = tmp_txt + '\n  ' + tmp_fp + ': ' 'NaN'

            subtext[tmp_sig] = tmp_txt

        tmp_error = [value for value in d_error.values() if not np.isnan(value)]

        subtext['mae'] = np.round(np.mean(np.abs(tmp_error)), 2)
        subtext['std'] = np.round(np.nanstd(tmp_error), 2)

        canvas = plot_fiducials(s=s, fp=fp1, savefig=False, show_fig=False, print_flag=False, use_tk=False, new_fig=True,
                                marker=marker1, title='', legend_loc='upper right', legend_fontsize=legend_fontsize,
                                marker_size=90, facecolor=False, subtext=subtext)

        if compare:
            title = 'Ref.1 & Ref.2'
            savingfolder = 'temp_dir'+os.sep+dname+os.sep+annot1+'-'+annot2
        else:
            title = 'Ref.  &  Det.'
            savingfolder = 'temp_dir'+os.sep+dname+os.sep+annot1+'-pyPPG'


        plot_fiducials(s=s, fp=fp2, savefig=True, show_fig=False, print_flag=False, use_tk=False, new_fig=False,
                       marker=marker2, title=title, legend_loc='upper right', legend_fontsize=legend_fontsize,
                       marker_size=50, facecolor=True, subtext=subtext, canvas=canvas, savingfolder=savingfolder)

        if not(show_fig): plt.close()

    ###########################################################################
    ########################### Print Distance Error ##########################
    ###########################################################################
    def print_error(self, d_error):
        # Print distance error
        print('\n', self.name, ' DIFF:')

        # Replace float values with integers in the dictionary
        for key, value in d_error.items():
            if isinstance(value, float):
                d_error[key] = int(value) if not np.isnan(value) else 'NaN'

        # Create a DataFrame from the modified dictionary
        df = pd.DataFrame(list(d_error.items()), columns=['Key', 'Value'])

        # Transpose the DataFrame
        transposed_df = df.transpose()

        # Print only the keys and values without column index
        print(transposed_df.to_string(index=False, header=False))

    ###########################################################################
    ################################# Save Data ###############################
    ###########################################################################
    def save_all_data(self, save, compare, dist_error,annot_error,annot_error2,dname, annot1, annot2):
        # Save distance error and annotation error in .mat file
        self.OutData[annot1+'_ae'] = annot_error           # 1st annotation error (0,1)
        self.OutData[annot1+'_fps'] = self.M_FID_1          # fiducial points of the 1st annotation
        if save:
            if compare:
                self.OutData[annot2+'_ae'] = annot_error2   # 2nd annotation error (0,1)
                self.OutData[annot2+'_fps'] = self.M_FID_2   # fiducial points of the 2nd annotation
                self.OutData[annot1 + '_'+annot2+'_diff'] = dist_error # annotation difference (sample)
                file_name = 'results'+os.sep+dname+os.sep+annot1+'-'+annot2+'.mat'
            else:
                self.OutData[annot1+'_pyPPG_diff'] = dist_error      # detection difference (sample)
                self.OutData['pyPPG_fps'] = self.M_FID_2   # fiducial points of the pyPPG
                try:
                    if len(self.annot_path2) > 0:
                        file_name = 'results'+os.sep+dname+os.sep+annot1+'-pyPPG'+'.mat'
                except:
                    file_name = 'results'+os.sep+dname+os.sep+annot2+'-pyPPG'+'.mat'

            # Convert DataFrames to MATLAB tables with headers
            matlab_data = {key: value.to_records(index=False) for key, value in self.OutData.items()}

            # Save the dictionary of MATLAB tables to a .mat file
            scipy.io.savemat(file_name, matlab_data, format='5', appendmat=False)
        return

    ###########################################################################
    ############################ Get PPG Validation ###########################
    ###########################################################################

    def get_validation(self, s=pyPPG.PPG, ref_fp=pd.DataFrame(), plt_sig=True, correction=pd.DataFrame(),dname='',annot1='', annot2=''):
        ## Create a PPG class
        s.correct = True
        s_class = PPG(s, check_ppg_len=False)

        ## Create a fiducial class
        fpex = FP.FpCollection(s_class)

        # Extract fiducial points
        peak = [ref_fp.sp[0]]
        onsets = [ref_fp.on[0], ref_fp.off[0]]

        ppg_fp = pd.DataFrame()
        ppg_fp['dn'] = np.array(fpex.get_dicrotic_notch(peak, onsets))

        det_dn = np.array(fpex.get_dicrotic_notch(peak, onsets))
        vpg_fp = fpex.get_vpg_fiducials(onsets)
        apg_fp = fpex.get_apg_fiducials(onsets, peak)
        jpg_fp = fpex.get_jpg_fiducials(onsets, apg_fp)

        ppg_fp['dp'] = fpex.get_diastolic_peak(onsets, ppg_fp.dn, apg_fp.e)

        if apg_fp.a[0] > 75:
            win_on = 75
        else:
            try:
                win_on = int(apg_fp.a[0])
            except:
                win_on = []
                pass

        try:
            ppg_fp['on'] = np.argmax(s.jpg[int(apg_fp.a[0]) - win_on:int(apg_fp.a[0])]) + apg_fp.a[0] - win_on
        except:
            ppg_fp['on'] = np.nan

        ppg_fp['off'] = ref_fp.off[0]
        ppg_fp['sp'] = np.argmax(s.ppg[int(ref_fp.on[0]):int(ref_fp.off[0])]) + ref_fp.on[0]

        # Merge fiducials
        det_fp = self.merge_fiducials(ppg_fp, vpg_fp, apg_fp, jpg_fp)

        # Correct fiducials
        det_fp = fpex.correct_fiducials(fiducials=det_fp, correction=correction)

        # Calculate distance error
        compare=False
        d_error = self.get_dist_error(ref_fp, det_fp, compare)

        # Plot fiducial points
        if plt_sig: self.plot_pulse_wave(s=s, fp1=ref_fp, fp2=det_fp, d_error=d_error, compare=compare, show_fig=False,dname=dname,annot1=annot1, annot2=annot2)

        return det_fp, d_error

    ###########################################################################
    ######################## Get Annotation Difference  #######################
    ###########################################################################
    def get_annot_diff(self, s, ref_fp,annot_path2,name,fid_names,compare,plt_sig,dname,annot1, annot2):
        ref1_fp = ref_fp

        annot_file2 = annot_path2 + os.sep + name + '.mat'
        annotf2 = scipy.io.loadmat(annot_file2)

        ref2_fp, annot_err = self.get_ref_fp(name, s.fs, fid_names, annotf2)
        ref2_fp = pd.DataFrame(ref2_fp, index=[0])

        # Calculate distance error
        annot_diff = self.get_dist_error(ref1_fp, ref2_fp, compare)

        # Plot fiducial points
        if plt_sig: self.plot_pulse_wave(s=s, fp1=ref1_fp, fp2=ref2_fp, d_error=annot_diff, compare=compare, show_fig=False,dname=dname,annot1=annot1, annot2=annot2)

        return ref2_fp, annot_diff, annot_err


    ###########################################################################
    ########################## Run PPG-BP Evaluation  #########################
    ###########################################################################

    def run_ppg_bp_eval(self, compare=False, plt_sig=False, save=False, prnt_e=True, correction=pd.DataFrame(), annot1='', annot2='', version='06', dname='YYYY_MM_DD_HH_MM'):
        # Define input directories and files
        ppg_sig_dir = 'PPG-BP_annot'
        ppg_file = os.sep+'PPG-BP_ref1.mat'
        annot_path = ppg_sig_dir+os.sep+annot1+version+'_PPG-BP_annot'+os.sep+'merged'
        self.annot_path2 = ppg_sig_dir + os.sep+annot2+version+'_PPG-BP_annot'+os.sep+'merged'
        sig_path = ppg_sig_dir + ppg_file
        input_sig = scipy.io.loadmat(sig_path)

        # Define output variables
        self.OutData = {}
        set_len = input_sig['ppg_data'].size

        fid_names = ('sp', 'on', 'dn', 'dp', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')

        if compare:
            f_names = fid_names[:]
        else:
            f_names = fid_names[:]

        dist_error = pd.DataFrame(columns=f_names)
        self.M_FID_1 = pd.DataFrame(columns=f_names)
        self.M_FID_2 = pd.DataFrame(columns=f_names)
        for n in f_names:
            exec(n + "=[]")
            exec(n + "r=[]")
            exec("dist_" + n + "=[]")
            temp_v = np.empty(set_len)
            temp_v[:] = np.NaN
            dist_error[n] = temp_v
            self.M_FID_1[n] = temp_v
            self.M_FID_2[n] = temp_v

        annot_error = pd.DataFrame(columns=f_names)
        for n in fid_names:
            annot_error[n] = temp_v

        annot_error2 = copy.deepcopy(annot_error)
        ppg_names={}

        for i in range(0, set_len):
            # Define sampling frequency, load filtered signal, 1st and 2nd derivative
            fs = input_sig['ppg_data']['fs'][0, 0][0][0]
            fs = np.squeeze(fs)

            # Load signal
            ppg_v = input_sig['ppg_data']['sig'][0, i]
            ppg_v = np.squeeze(ppg_v)

            # Load annotated fiducial points
            self.name = input_sig['ppg_data']['name'][0, i][0]
            ppg_names[i]=self.name
            annot_file = annot_path + os.sep + self.name + '.mat'
            annot = scipy.io.loadmat(annot_file)

            # Get reference annotation
            ref_fp, error = self.get_ref_fp(self.name, fs, fid_names, annot)
            self.M_FID_1.iloc[i] = ref_fp
            ref_fp = pd.DataFrame(ref_fp, index=[0])
            annot_error.iloc[i] = error

            # Create struct for PPG signal
            s = DotMap()
            s.fs = fs
            s.name = self.name
            s.v = ppg_v
            s.ppg = ppg_v

            ## Preprocessing
            prep = PP.Preprocess()

            # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
            s.filtering = True
            s.ppg, s.vpg, s.apg, s.jpg = prep.get_signals(s)


            if compare:
                ref2_fp, annot_diff, annot_err = self.get_annot_diff(s, ref_fp, self.annot_path2, self.name, fid_names, compare, plt_sig,dname,annot1, annot2)
                dist_error.iloc[i] = annot_diff
                annot_error2.iloc[i] = annot_err
                self.M_FID_2.iloc[i] = ref2_fp[list(f_names)]
                self.save_all_data(save, compare, dist_error, annot_error, annot_error2,dname, annot1,annot2)
            else:
                det_fp, annot_diff = self.get_validation(s=s, ref_fp=ref_fp, plt_sig=plt_sig, correction=correction,dname=dname,annot1=annot1, annot2=annot2)
                dist_error.iloc[i] = annot_diff
                self.M_FID_2.iloc[i] = det_fp[list(f_names)]
                self.save_all_data(save, compare, dist_error, annot_error,annot_error2, dname, annot1, annot2)

            # Print distance error
            if prnt_e: self.print_error(annot_diff)

        # Calculate Mean Absolute Error (MAE), Standard Deviation (STD), and Mean Error (BIAS)
        MAE = {}
        STD = {}
        BIAS = {}
        for n in f_names:
            exec("MAE['" + n + "'] = np.round(np.nanmean(np.absolute(dist_error." + n + ")),2)")
            exec("STD['" + n + "'] = np.round(np.nanstd(dist_error." + n + "), 2)")
            exec("BIAS['" + n + "'] = np.round(np.nanmean(dist_error." + n + "),2)")

        # Insert statistics into the dataframe
        dist_error.index = dist_error.index + 3
        dist_error.loc[0] = MAE
        dist_error.loc[1] = STD
        dist_error.loc[2] = BIAS
        dist_error=dist_error.sort_index()
        dist_error = dist_error.rename({0: 'MAE'})
        dist_error = dist_error.rename({1: 'STD'})
        dist_error = dist_error.rename({2: 'BIAS'})
        dist_error.index=dist_error.index.map(lambda x: ppg_names[x-3] if isinstance(x, int) else x)

        # Save Results
        date = datetime.now()
        if compare:
            fname=annot1+'-'+annot2
        else:
            fname=annot1+'-pyPPG'

        path1='results'+os.sep+dname+os.sep+fname+'_diffs.csv'
        dist_error.to_csv(path1, index=True)

        path2 = 'results' + os.sep + dname + os.sep + fname + '_params.txt'
        params='fH: '+str(prep.fH)+'\nfL: '+str(prep.fL)+'\norder: '+str(prep.order)+'\nsm_win: '+str(prep.sm_wins)
        # Open the file in write mode ('w') and write the string
        with open(path2, 'w') as file:
            file.write(params)

        # Print results
        print('-------------------------------------------')
        print('MAE', ': ', MAE)
        print('STD', ': ', STD)
        print('BIAS', ': ', BIAS)
        print('Program finished!')

    ###########################################################################
    ####################### Extract Pulse Wave Features  ######################
    ###########################################################################
    def pw_extraction(self, data_path="", filtering=True, fL=0, fH=12, order=4, sm_wins={'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10},
                       correction=pd.DataFrame(), savefig=True, savingfolder='', show_fig=False, print_flag=True):

        ## Loading a raw PPG signal
        signal = load_data(data_path=data_path)
        # signal.fs=300

        ## Preprocessing
        # Initialise the filters
        prep = PP.Preprocess(fL=fL, fH=fH, order=order, sm_wins=sm_wins)

        # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
        signal.filtering = filtering
        signal.ppg, signal.vpg, signal.apg, signal.jpg = prep.get_signals(s=signal)

        ## Create a PPG class
        signal.correction = correction
        s = PPG(signal, check_ppg_len=False)

        ## Create a fiducial class
        fpex = FP.FpCollection(s)

        # Extract fiducial points
        ppg_fp = pd.DataFrame()
        ppg_fp['on'] = [0]
        ppg_fp['off'] = [len(signal.ppg)-1]
        ppg_fp['sp'] = [np.argmax(signal.ppg)]

        peak = [ppg_fp.sp.iloc[0]]
        onsets = [ppg_fp.on.iloc[0], ppg_fp.off.iloc[0]]

        ppg_fp['dn'] = np.array(fpex.get_dicrotic_notch(peak, onsets))

        det_dn = np.array(fpex.get_dicrotic_notch(peak, onsets))
        vpg_fp = fpex.get_vpg_fiducials(onsets)
        apg_fp = fpex.get_apg_fiducials(onsets, peak)
        jpg_fp = fpex.get_jpg_fiducials(onsets, apg_fp)

        ppg_fp['dp'] = fpex.get_diastolic_peak(onsets, ppg_fp.dn, apg_fp.e)

        det_fp = self.merge_fiducials(ppg_fp, vpg_fp, apg_fp, jpg_fp)

        if bool(len(correction)): det_fp = fpex.correct_fiducials(fiducials=det_fp, correction=correction)

        det_fp = det_fp.fillna(0).astype(int)

        det_fp_new = pd.concat([det_fp,pd.DataFrame({'on': [ppg_fp.off.iloc[0]]})])

        # Replace float values with integers in the dictionary
        det_fp_new = det_fp_new.fillna(0).astype(int)
        det_fp_new.index = [0,1]

        if print_flag:print(det_fp)

        ## Plot fiducial points
        plot_fiducials(s=s, fp=det_fp, savefig=savefig, savingfolder=savingfolder, show_fig=show_fig,
                       print_flag=print_flag, legend_fontsize=9, marker_size=90, use_tk=False)

        bm = self.get_pw_bm(s=s, fp=det_fp_new)

        return s, det_fp, bm

    ###########################################################################
    ###################### Extract Pulse Wave Biomarkers  #####################
    ###########################################################################
    def get_pw_bm(self, s=pyPPG.PPG, fp=pd.DataFrame()):
        ## Get Biomarkers and Statistics
        # Initialise the biomarkers package
        fp = Fiducials(fp=fp)
        bmex = BM.BmCollection(s=s, fp=fp)

        # Extract biomarkers
        bm_defs, bm_vals = bmex.get_biomarkers(get_stat=False)

        # Create a biomarkers class
        bm = Biomarkers(bm_defs=bm_defs, bm_vals=bm_vals)

        return bm

    ###########################################################################
    ###################### Extract Pulse Wave Biomarkers  #####################
    ###########################################################################
    def extact_pw_feat(self, number_of_rec):

        savingfolder = 'temp_dir' + os.sep + 'PW_anal_01'
        all_fp = pd.DataFrame()
        all_bm_vals = {}

        for i in range(0, number_of_rec):
            name = "single_pw_sample"

            data_path = 'Single_PW' + os.sep + name + '.mat'
            s, fp, bm = pwex.pw_extraction(data_path=data_path, filtering=True, fL=0, fH=12, order=4,
                                           sm_wins={'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10}, correction=correction,
                                           savefig=True, savingfolder=savingfolder, show_fig=False, print_flag=True)

            all_fp = pd.concat([all_fp, fp])

            for bm_key in bm.bm_vals.keys():
                if i == 0:
                    all_bm_vals[bm_key] = bm.bm_vals[bm_key]
                else:
                    all_bm_vals[bm_key] = pd.concat([all_bm_vals[bm_key], bm.bm_vals[bm_key]])

        for bm_key in bm.bm_vals.keys():
            all_bm_vals[bm_key].index = range(0, number_of_rec)

        all_fp.index = range(0, number_of_rec)

        bm_stats = get_statistics(all_fp.sp, all_fp.on, all_bm_vals)

        BMs = Biomarkers(bm_defs=bm.bm_defs, bm_vals=all_bm_vals, bm_stats=bm_stats)
        FPs = Fiducials(fp=all_fp)

        ## Save data
        save_data(s=s, fp=FPs, bm=BMs, savingformat='csv', savingfolder=savingfolder, print_flag=True)

###########################################################################
############################ MAIN PPG Analysis  ###########################
###########################################################################
if __name__ == '__main__':

    # Flag for package usage
    ppgbp=True      # validation of PPG-BP dataset
    pw_ext=False    # extract features of pulse waves
    plts=True      # plot signal

    # Initialise the pulse wave package
    pwex = PulseWaveAnal()

    # Initialise the correction
    correction = pd.DataFrame()
    corr_on = ['on', 'dp', 'v', 'w', 'f']
    corr_off = ['dn']

    correction.loc[0, corr_on] = True
    correction.loc[0, corr_off] = False

    date = datetime.now()
    dname = str(date.year) + '_' + str(date.month) + '_' + str(date.day) + '_' + str(date.hour) + '_' + str(date.minute)
    os.makedirs('results' + os.sep + dname, exist_ok=True)

    # Run PPG-BP Evaluation
    if  ppgbp:
        pwex.run_ppg_bp_eval(compare=False, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='MG', annot2='PC', version='06', dname=dname)
        pwex.run_ppg_bp_eval(compare=False, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='PC', annot2='MG', version='06', dname=dname)
        pwex.run_ppg_bp_eval(compare=True, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='MG', annot2='PC', version='06', dname=dname)

    # Extract Pulse Wave Features
    if pw_ext:
        pwex.extact_pw_feat(number_of_rec=10)

    print('End of analysis!')




