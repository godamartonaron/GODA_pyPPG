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
import pkg_resources

import json
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
                exec("ref_fp['on'] = ref_on[0]")
            elif n=="off":
                exec("ref_fp['off'] = ref_on[1]")
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

    def plot_pulse_wave(self, s=pd.DataFrame(),fp1=pd.DataFrame(),fp2=pd.DataFrame(),d_error={},compare=False,show_fig=False,dname='',annot1='', annot2='',detector=''):
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
            title = annot1+' & '+annot2
            savingfolder = 'results'+os.sep+dname+os.sep+annot1+'_'+annot2
        else:
            title = annot1+' & '+detector
            tmp_dir='results' + os.sep + dname + os.sep + detector+ os.sep
            os.makedirs(tmp_dir, exist_ok=True)
            savingfolder = tmp_dir+annot1+'_'+detector


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
        #self.OutData[annot1+'_ae'] = annot_error           # 1st annotation error (0,1)
        self.OutData[annot1+'_fps'] = self.M_FID_1          # fiducial points of the 1st annotation
        if save:
            if compare:
                tmp_dir = 'results'+os.sep+dname+os.sep+annot1+'_'+annot2+os.sep
                os.makedirs(tmp_dir, exist_ok=True)

                #self.OutData[annot2+'_ae'] = annot_error2   # 2nd annotation error (0,1)
                self.OutData[annot2+'_fps'] = self.M_FID_2   # fiducial points of the 2nd annotation
                self.OutData[annot1 + '_'+annot2+'_diff'] = dist_error # annotation difference (sample)
                mat_file = tmp_dir+annot1+'_'+annot2+'.mat'

            else:
                self.OutData[annot1+'_pyPPG_diff'] = dist_error      # detection difference (sample)
                self.OutData['pyPPG_fps'] = self.M_FID_2   # fiducial points of the pyPPG

                tmp_dir1 = 'results' + os.sep + dname + os.sep + 'pyPPG'+ os.sep+annot1+'_pyPPG' + os.sep
                os.makedirs(tmp_dir1, exist_ok=True)

                tmp_dir2 = 'results' + os.sep + dname + os.sep + 'pyPPG'+ os.sep+annot2+'_pyPPG' + os.sep
                os.makedirs(tmp_dir2, exist_ok=True)

                try:
                    if len(self.annot_path2) > 0:
                        mat_file = tmp_dir1 + annot1 + '_pyPPG'+'.mat'
                except:
                    mat_file = tmp_dir2 + annot2+'_pyPPG'+'.mat'

            # Convert DataFrames to MATLAB tables with headers
            matlab_data = {key: value.to_records(index=False) for key, value in self.OutData.items()}

            # Save the dictionary of MATLAB tables to a .mat file
            scipy.io.savemat(mat_file, matlab_data, format='5', appendmat=False)
        return

    ###########################################################################
    ############################ Get PPG Validation ###########################
    ###########################################################################

    def get_validation(self, s=pyPPG.PPG, ref_fp=pd.DataFrame(), plt_sig=True, correction=pd.DataFrame(),dname='',annot1='', annot2='', detector=''):
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
        if plt_sig: self.plot_pulse_wave(s=s, fp1=ref_fp, fp2=det_fp, d_error=d_error, compare=compare, show_fig=False,dname=dname,annot1=annot1, annot2=annot2, detector=detector)

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
    ############################ Get PPG statistics  ##########################
    ###########################################################################
    def get_stats(self,fp_names,dist_error,ppg_IDs, dname, params, filename, detector_dir, prnt):
        # Calculate Mean Absolute Error (MAE), Standard Deviation (SD), and Mean Error (BIAS)
        MAE = {}
        SD = {}
        BIAS = {}
        for n in fp_names:
            exec("MAE['" + n + "'] = np.round(np.nanmean(np.absolute(dist_error." + n + ")),2)")
            exec("SD['" + n + "'] = np.round(np.nanstd(dist_error." + n + "), 2)")
            exec("BIAS['" + n + "'] = np.round(np.nanmean(dist_error." + n + "),2)")

        # Insert statistics into the dataframe
        dist_error.index = dist_error.index + 3
        dist_error.loc[0] = MAE
        dist_error.loc[1] = SD
        dist_error.loc[2] = BIAS
        dist_error=dist_error.sort_index()
        dist_error = dist_error.rename({0: 'MAE'})
        dist_error = dist_error.rename({1: 'SD'})
        dist_error = dist_error.rename({2: 'BIAS'})
        dist_error.index=dist_error.index.map(lambda x: ppg_IDs[x-3] if isinstance(x, int) else x)

        tmp_dir = 'results'+os.sep+dname+os.sep+detector_dir+filename+ os.sep
        os.makedirs(tmp_dir, exist_ok=True)

        path=tmp_dir+filename+'_diffs.csv'
        dist_error.to_csv(path, index=True)

        if len(params)>0:
            path2 = tmp_dir+ filename + '_params.txt'
            # Open the file in write mode ('w') and write the string
            with open(path2, 'w') as file:
                file.write(params)


        if prnt:
            # Print results
            print('-------------------------------------------')
            print('MAE', ': ', MAE)
            print('SD', ': ', SD)
            print('BIAS', ': ', BIAS)

    ###########################################################################
    ########################## Run PPG-BP Evaluation  #########################
    ###########################################################################

    def run_ppg_bp_eval(self, compare=False, plt_sig=False, save=False, prnt_e=True, correction=pd.DataFrame(), annot1='', annot2='', version='VV', dname='YYYY_MM_DD_HH_MM', prnt=False):
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

        fp_names = ('sp', 'on', 'dn', 'dp', 'off', 'u', 'v', 'w', 'a', 'b', 'c', 'd', 'e', 'f', 'p1', 'p2')

        dist_error = pd.DataFrame(columns=fp_names)
        ref1_fps = pd.DataFrame(columns=fp_names)
        ref2_fps = pd.DataFrame(columns=fp_names)
        det_fps = pd.DataFrame(columns=fp_names)
        self.M_FID_1 = pd.DataFrame(columns=fp_names)
        self.M_FID_2 = pd.DataFrame(columns=fp_names)
        annot_error = pd.DataFrame(columns=fp_names)
        for n in fp_names:
            exec(n + "=[]")
            exec(n + "r=[]")
            exec("dist_" + n + "=[]")
            temp_v = np.empty(set_len)
            temp_v[:] = np.NaN
            dist_error[n] = temp_v
            ref1_fps[n] = temp_v
            ref2_fps[n] = temp_v
            det_fps[n] = temp_v
            self.M_FID_1[n] = temp_v
            self.M_FID_2[n] = temp_v
            annot_error[n] = temp_v

        annot_error2 = copy.deepcopy(annot_error)
        ppg_IDs={}

        for i in range(0, set_len):
            # Define sampling frequency, load filtered signal, 1st and 2nd derivative
            fs = input_sig['ppg_data']['fs'][0, 0][0][0]
            fs = np.squeeze(fs)

            # Load signal
            ppg_v = input_sig['ppg_data']['sig'][0, i]
            ppg_v = np.squeeze(ppg_v)

            # Load annotated fiducial points
            self.name = input_sig['ppg_data']['name'][0, i][0]
            ppg_IDs[i]=self.name
            annot_file = annot_path + os.sep + self.name + '.mat'
            annot = scipy.io.loadmat(annot_file)

            # Get reference annotation
            ref_fp, error = self.get_ref_fp(self.name, fs, fp_names, annot)
            self.M_FID_1.iloc[i] = ref_fp
            ref_fp = pd.DataFrame(ref_fp, index=[0])
            ref1_fps.iloc[i] = ref_fp
            annot_error.iloc[i] = error

            # Create struct for PPG signal
            s = DotMap()
            s.fs = fs
            s.name = self.name
            s.v = ppg_v
            s.ppg = ppg_v

            ## Preprocessing
            # sm_wins = {'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10}
            # prep = PP.Preprocess(fL=0, fH=12, order=4, sm_wins=sm_wins)
            prep = PP.Preprocess()

            # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
            s.filtering = True
            s.ppg, s.vpg, s.apg, s.jpg = prep.get_signals(s)

            sigs_dir='results' + os.sep + dname + os.sep + 'pyPPG' + os.sep + 'Signals' +os.sep
            os.makedirs(sigs_dir, exist_ok=True)
            tmp_mat_name=self.name+'.mat'
            path=sigs_dir+tmp_mat_name
            if not(os.path.exists(path)):
                scipy.io.savemat(path, s)

            if compare:
                ref2_fp, annot_diff, annot_err = self.get_annot_diff(s, ref_fp, self.annot_path2, self.name, fp_names, compare, plt_sig,dname,annot1, annot2)
                ref2_fps.iloc[i] = ref2_fp
                dist_error.iloc[i] = annot_diff
                annot_error2.iloc[i] = annot_err
                self.M_FID_2.iloc[i] = ref2_fp[list(fp_names)]
                self.save_all_data(save, compare, dist_error, annot_error, annot_error2,dname, annot1,annot2)
                detector_dir = ''
            else:
                det_fp, annot_diff = self.get_validation(s=s, ref_fp=ref_fp, plt_sig=plt_sig, correction=correction,dname=dname,annot1=annot1, annot2=annot2, detector='pyPPG')
                det_fps.iloc[i] = det_fp[list(fp_names)]
                dist_error.iloc[i] = annot_diff
                self.M_FID_2.iloc[i] = det_fp[list(fp_names)]
                self.save_all_data(save, compare, dist_error, annot_error,annot_error2, dname, annot1, annot2)
                detector_dir='pyPPG'+os.sep

            # Print distance error
            if prnt_e: self.print_error(annot_diff)

        if compare:
            filename=annot1+'_'+annot2
            tmp_dir = 'results' + os.sep + dname + os.sep + 'MG_PC' + os.sep
            csv_file1 = tmp_dir + annot1 + '_fps.csv'
            csv_file2 = tmp_dir + annot2 + '_fps.csv'
            ref1_fps = ref1_fps.set_index(pd.Index(ppg_IDs.values()))
            ref2_fps = ref2_fps.set_index(pd.Index(ppg_IDs.values()))
            ref1_fps.to_csv(csv_file1, index=True)
            ref2_fps.to_csv(csv_file2, index=True)
            corr_txt=''
        else:
            filename=annot1+'_pyPPG'
            tmp_dir = 'results' + os.sep + dname + os.sep + 'pyPPG' + os.sep
            csv_file1 = tmp_dir + annot1 + '_pyPPG' + os.sep + annot1 + '_fps.csv'
            csv_file2 = tmp_dir + annot1 + '_pyPPG' + os.sep + 'pyPPG_fps.csv'
            ref1_fps = ref1_fps.set_index(pd.Index(ppg_IDs.values()))
            det_fps = det_fps.set_index(pd.Index(ppg_IDs.values()))
            ref1_fps = ref1_fps.rename_axis('ID')
            det_fps = det_fps.rename_axis('ID')
            ref1_fps.to_csv(csv_file1, index=True)
            det_fps.to_csv(csv_file2, index=True)
            corr_txt = '\ncorrection: '+ str(correction.to_dict()).replace("{0:", "").replace("},", ",").replace("}}", "}")

        # Get the version using pkg_resources
        package_name = 'pyPPG'
        pyPPG_version = pkg_resources.get_distribution(package_name).version
        params = 'pyPPG version: ' + pyPPG_version + '\nfH: ' + str(prep.fH) + ' Hz\nfL: ' + str(prep.fL) + ' Hz\norder: ' + str(prep.order) + '\nsm_win: ' + str(prep.sm_wins)  + corr_txt

        path = 'results' + os.sep + dname + os.sep + 'ppg_IDs.csv'
        df_ppg_IDs= pd.DataFrame.from_dict(ppg_IDs, orient='index', columns=['ID'])
        df_ppg_IDs=df_ppg_IDs.set_index(df_ppg_IDs.index + 1)
        df_ppg_IDs.to_csv(path, index=True)

        self.get_stats(fp_names, dist_error, ppg_IDs, dname, params, filename, detector_dir, prnt)

    ###########################################################################
    ##################### Benchmark PPG-BP fiducial points  ###################
    ###########################################################################

    def benchmark_PPG_BP(self, detector, dname, plt, prnt):
        # set loading directory
        tmp_dir='results' + os.sep + dname + os.sep+ detector +os.sep

        for annotator in ['MG','PC']:
            name=annotator+'_'+detector
            input_filename = tmp_dir + name+os.sep+name + '.mat'
            tmp_data=scipy.io.loadmat(input_filename)

            # get distance error
            diff_name = name + '_diff'
            raw_diffs = tmp_data[diff_name]
            flat_diffs  = raw_diffs.flatten()
            dist_error = pd.DataFrame(flat_diffs.tolist(), columns=raw_diffs.dtype.names)
            dist_error = dist_error.applymap(lambda x: np.squeeze(x))

            # get detection
            det_name = detector + '_fps'
            raw_dets = tmp_data[det_name]
            flat_dets  = raw_dets.flatten()
            detection = pd.DataFrame(flat_dets.tolist(), columns=raw_dets.dtype.names)
            detection = detection.applymap(lambda x: np.squeeze(x))

            # get annotation
            ref_name = annotator + '_fps'
            raw_refs = tmp_data[ref_name]
            flat_refs  = raw_refs.flatten()
            reference = pd.DataFrame(flat_refs.tolist(), columns=raw_refs.dtype.names)
            reference = reference.applymap(lambda x: np.squeeze(x))

            # load PPG IDs
            fp_names=raw_diffs.dtype.names
            path = 'results' + os.sep + dname + os.sep + 'ppg_IDs.csv'
            ppg_IDs=pd.read_csv(path).ID

            if detector=="PPGFeat":
                params = 'fH: ' + '12 Hz' + '\nfL: ' + '0.5 Hz'
            elif detector=="PulseAnal" :
                sm_wins = {'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10}
                prep = PP.Preprocess(fL=0, fH=12, order=4, sm_wins=sm_wins)
                params = 'fH: ' + str(prep.fH) + ' Hz\nfL: ' + str(prep.fL) + ' Hz\norder: ' + str(prep.order) + '\nsm_win: ' + str(prep.sm_wins)
            else:
                params = ''

            if plt:
                sigs_dir = 'results' + os.sep + dname + os.sep + 'pyPPG' + os.sep + 'Signals' + os.sep

                n = 0
                for ID in ppg_IDs:
                    print(annotator + '-' +detector+' '+ str(n + 1) + '/219' + ': ' + ID)
                    # Plot fiducial points
                    tmp_mat_name = ID + '.mat'
                    path = sigs_dir + tmp_mat_name
                    sc = scipy.io.loadmat(path)
                    keys = list(sc.keys())[3:]
                    s = DotMap((key, sc[key][0]) for key in keys)
                    plt_sig = 1
                    r = dict(reference.iloc[n].apply(pd.to_numeric, errors='coerce').astype('Int32'))
                    ref_fp = pd.DataFrame([{key: np.NaN if pd.isna(value) else value for key, value in r.items()}])
                    d = dict(detection.iloc[n].apply(pd.to_numeric, errors='coerce').astype('Int32'))
                    det_fp = pd.DataFrame([{key: np.NaN if pd.isna(value) else value for key, value in d.items()}])
                    if plt_sig: self.plot_pulse_wave(s=s, fp1=ref_fp, fp2=det_fp, d_error=dist_error.iloc[n].to_dict(),
                                                     compare=False, show_fig=False, dname=dname, annot1=annotator, annot2='',
                                                     detector=detector)
                    n = n + 1

            detector_dir=detector+os.sep
            self.get_stats(fp_names, dist_error, ppg_IDs, dname, params, name, detector_dir, prnt)

            tmp_dir = 'results' + os.sep + dname + os.sep + detector + os.sep
            csv_file1 = tmp_dir + annotator + '_' + detector + os.sep + annotator + '_fps.csv'
            csv_file2 = tmp_dir + annotator + '_' + detector + os.sep + detector+'_fps.csv'
            ref_fps = reference.set_index(pd.Index(ppg_IDs.values))
            det_fps = detection.set_index(pd.Index(ppg_IDs.values))
            ref_fps = ref_fps.rename_axis('ID')
            det_fps = det_fps.rename_axis('ID')
            ref_fps.to_csv(csv_file1, index=True)
            det_fps.to_csv(csv_file2, index=True)

        MG_diff_path = tmp_dir + 'MG_' + detector +os.sep + 'MG_' + detector + '_diffs.csv'
        PC_diff_path = tmp_dir + 'PC_' + detector +os.sep + 'PC_' + detector + '_diffs.csv'
        MG_diff=pd.read_csv(MG_diff_path, sep=',').head(3)
        PC_diff = pd.read_csv(PC_diff_path, sep=',').head(3)

        mean_diff=(MG_diff[list(fp_names)]+PC_diff[list(fp_names)])/2
        row_names = MG_diff['Unnamed: 0'].tolist()
        mean_diff.insert(0, 'Stats', row_names)

        print('---------------------------------')
        print(detector+' results:')
        print(mean_diff)

        path = tmp_dir+detector+'_results.csv'
        mean_diff.to_csv(path, index=False)
        print('Output path: ',path)

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
    def extract_pw_feat(self, datafolder, savingfolder,correction,savefig,show_fig):

        all_fp = pd.DataFrame()
        all_bm_vals = {}
        all_pw = os.listdir(datafolder)
        number_of_rec=len(all_pw)
        for i in range(0,number_of_rec):
            name = all_pw[i]

            data_path = datafolder + os.sep + name
            s, fp, bm = self.pw_extraction(data_path=data_path, filtering=True, fL=0, fH=12, order=4,
                                           sm_wins={'ppg': 50, 'vpg': 10, 'apg': 10, 'jpg': 10}, correction=correction,
                                           savefig=savefig, savingfolder=savingfolder, show_fig=show_fig, print_flag=True)

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
    ########################### Evaluation of PPG-BP  #########################
    ###########################################################################
    def eval_PPG_BP(self,plts,correction,dname,prnt):
        os.makedirs('results' + os.sep + dname, exist_ok=True)

        # Run PPG-BP Evaluation
        self.run_ppg_bp_eval(compare=False, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='MG',
                             annot2='PC', version='final', dname=dname, prnt=prnt)
        self.run_ppg_bp_eval(compare=False, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='PC',
                             annot2='MG', version='final', dname=dname, prnt=prnt)
        self.run_ppg_bp_eval(compare=True, plt_sig=plts, save=True, prnt_e=True, correction=correction, annot1='MG',
                             annot2='PC', version='final', dname=dname, prnt=prnt)

    ###########################################################################
    ######################## Extract Pulse Wave Features  #####################
    ###########################################################################
    def pw_extraction(self,correction):
        savingfolder = 'temp_dir' + os.sep + 'PW_anal_01'
        datafolder = 'Single_PW'
        self.extract_pw_feat(datafolder, savingfolder,correction,savefig=False,show_fig=False)

    ###########################################################################
    #############################  Run Benchmarking  ##########################
    ###########################################################################
    def run_benchmarking(self,detector, dname, plt, prnt):
        self.benchmark_PPG_BP(detector, dname, plt, prnt)

###########################################################################
############################ MAIN PPG Analysis  ###########################
###########################################################################
if __name__ == '__main__':
    # Initialise the pulse wave package
    pwex = PulseWaveAnal()

    # Initialise the correction
    correction = pd.DataFrame()
    corr_on = ['on', 'v', 'w', 'f']
    corr_off = ['dn', 'dp']
    correction.loc[0, corr_on] = True
    correction.loc[0, corr_off] = False

    # Save time for the results of the output folders
    date = datetime.now()
    dname = str(date.year) + '_' + str(date.month) + '_' + str(date.day) + '_' + str(date.hour) + '_' + str(date.minute)

    # Run PPG-BP Evaluation
    pwex.eval_PPG_BP(plts=False, correction=correction, dname=dname)

    # Run Benchmarking
    pwex.benchmark_PPG_BP(detector='PPGFeat', dname='2024_1_7_16_1',plt=False, prnt=False)
    pwex.benchmark_PPG_BP(detector='pyPPG', dname='2024_1_7_16_1',plt=False, prnt=False)
    pwex.benchmark_PPG_BP(detector='PulseAnal', dname='2024_1_7_16_1',plt=False, prnt=False)

    # Extract Pulse Wave Features
    #pwex.pw_extraction(correction=correction)

    print('End of analysis!')




