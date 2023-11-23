import pyPPG
from pyPPG import PPG, Fiducials, Biomarkers
from pyPPG.datahandling import load_data, plot_fiducials, save_data, load_fiducials
import pyPPG.preproc as PP
import pyPPG.fiducials as FP
import pyPPG.biomarkers as BM
import pyPPG.ppg_sqi as SQI

import numpy as np
import sys
import json
import pandas as pd
import scipy.io

###########################################################################
################################## EXAMPLE ################################
###########################################################################
def ppg_example(data_path="", fs=0, start_sig=0, end_sig=-1, fiducials=pd.DataFrame(), process_type="both", channel="Pleth",
                filtering=True, fL=0.5000001, fH=12, order=4, sm_wins={'ppg':50,'vpg':10,'apg':10,'jpg':10}, correction=pd.DataFrame(),
                savingfolder="temp_dir", savefig=True, show_fig=True, savingformat="both", print_flag=True, use_tk=False, check_ppg_len=True,
                saved_fiducials="", savedata=True):
    '''
    This is an example code for PPG analysis. The main parts:
        1) Loading a raw PPG signal: various file formats such as .mat, .csv, .txt, or .edf.
        2) Get Fiducial points: extract the fiducial points of PPG, PPG', PPG'' and PPG'" signals
        3) Plot Fiducial Points
        4) Get Biomarkers: extract 74 PPG biomarkers in four categories:
            - PPG signal
            - Signal ratios
            - PPG derivatives
            - Derivatives ratios
        5) Get Statistics: summary of the 74 PPG biomarkers
        6) SQI calculation: calculates the PPG Signal Quality Index
        7) Save data: save the extracted Fiducial points, Biomarkers, and Statistics into .csv file

    :param data_path: path of the PPG signal
    :type data_path: str
    :param fs: sampling_frequency
    :type fs: int
    :param start_sig: beginning the of signal in sample
    :type start_sig: int
    :param end_sig: end of the signal in sample
    :type end_sig: int
    :param fiducials: DataFrame of the fiducial points
    :type fiducials: pyPPG.Fiducials DataFrame
    :param process_type: the type of the process, which can be "fiducials", "biomarkers", or "both"
    :type process_type: str
    :param channel: channel of the .edf file
    :type channel: channel of the .edf file
    :param filtering: a bool for filtering
    :type filtering: bool
    :param fL: Lower cutoff frequency (Hz)
    :type fL: float
    :param fH: Upper cutoff frequency (Hz)
    :type fH: float
    :param order: Filter order
    :type order: int
    :param sm_wins: dictionary of smoothing windows in millisecond:
        - ppg: window for PPG signal
        - vpg: window for PPG' signal
        - apg: window for PPG" signal
        - jpg: window for PPG'" signal
    :type sm_wins: dict
    :param correction: DataFrame where the key is the name of the fiducial points and the value is bool
    :type correction: DataFrame
    :param savingfolder: location of the saved data
    :type savingfolder: str
    :param savefig: a bool for current figure saving
    :type savefig: bool
    :param show_fig: a bool for show figure
    :type show_fig: bool
    :param savingformat: file format of the saved date, the provided file formats .mat, .csv, or both
    :type savingformat: str
    :param print_flag: a bool for print message
    :type print_flag: bool
    :param use_tk: a bool for using tkinter interface
    :type use_tk: bool
    :param check_ppg: a bool for checking ppg length and sampling frequency
    :type check_ppg: bool
    :param saved_fiducials: path of the file of the saved fiducial points
    :type saved_fiducials: str
    :param savedata: a bool for saving data
    :type savedata: bool

    :return: file_names: dictionary of the saved file names

    Example:

        .. code-block:: python

            from pyPPG.example import ppg_example

            # run example code
            ppg_example()

    '''

    ## Loading a raw PPG signal
    signal = load_data(data_path=data_path, fs=fs, start_sig=start_sig, end_sig=end_sig, channel=channel, use_tk=True, print_flag=print_flag)

    ## Preprocessing
    # Initialise the filters
    prep = PP.Preprocess(fL=fL, fH=fH, order=order, sm_wins=sm_wins)

    # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
    signal.filtering = filtering
    signal.fL = fL
    signal.fH = fH
    signal.order = order
    signal.sm_wins = sm_wins
    signal.ppg, signal.vpg, signal.apg, signal.jpg = prep.get_signals(s=signal)

    # Initialise the correction for fiducial points
    corr_on = ['on', 'dn', 'dp', 'v', 'w', 'f']
    correction.loc[0, corr_on] = True
    signal.correction=correction

    ## Create a PPG class
    s = PPG(s=signal, check_ppg_len=check_ppg_len)

    ## Get Fiducial points
    if process_type == 'fiducials' or process_type == 'both':
        # Initialise the fiducials package
        fpex = FP.FpCollection(s=s)

        # Extract fiducial points
        fiducials = fpex.get_fiducials(s=s)
        if print_flag: print("Fiducial points:\n", fiducials + s.start_sig)

        # Create a fiducials class
        fp = Fiducials(fp=fiducials)

        # Save data
        if savedata:
            fp_new = Fiducials(fp=fp.get_fp() + s.start_sig)
            file_names=save_data(savingformat=savingformat, savingfolder=savingfolder, print_flag=print_flag, s=s, fp=fp_new)

    ## PPG SQI

        # Calculate SQI
        ppgSQI = round(np.mean(SQI.get_ppgSQI(ppg=s.ppg, fs=s.fs, annotation=fp.sp)) * 100, 2)
        if print_flag: print('Mean PPG SQI: ', ppgSQI, '%')

    ## Plot fiducial points
        plot_fiducials(s=s, fp=fp, savefig=savefig, savingfolder=savingfolder, show_fig=show_fig, print_flag=print_flag, use_tk=use_tk)

    ## Load saved fiducial points from MATLAB struct
    if ".mat" in saved_fiducials:
        tmp_fp1 = load_fiducials(saved_fiducials=saved_fiducials)
        tmp_fp2 = tmp_fp1[(tmp_fp1['on']>= s.start_sig) & (tmp_fp1['off']<= s.end_sig)]
        fiducials = tmp_fp2-s.start_sig
        fiducials.index =range(0,len(fiducials))

        ## Get Biomarkers and Statistics
    if (process_type == 'biomarkers' or process_type == 'both') and len(fiducials)>0:
        # Initialise the biomarkers package
        fp = Fiducials(fp=fiducials)

        bmex = BM.BmCollection(s=s, fp=fp)

        # Extract biomarkers
        bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()

        if print_flag:
            tmp_keys = bm_stats.keys()
            print('Statistics of the biomarkers:')
            for i in tmp_keys: print(i, '\n', bm_stats[i])

        # Create a biomarkers class
        bm = Biomarkers(bm_defs=bm_defs, bm_vals=bm_vals, bm_stats=bm_stats)

        # Save data
        if savedata:
            fp_new = Fiducials(fp=fp.get_fp() + s.start_sig)
            file_names=save_data(savingformat=savingformat, savingfolder=savingfolder, print_flag=print_flag, s=s, fp=fp_new, bm=bm)

    if print_flag: print('Program finished')

    return file_names


###########################################################################
############################## RUN EXAMPLE CODE ###########################
###########################################################################
if __name__ == "__main__":

    if len(sys.argv) > 1:
        input_data = json.loads(sys.argv[1])
        function_name = input_data['function']
        function_args = input_data['args']

        if function_name == 'ppg_example':
            file_names = ppg_example(**function_args)
            print(json.dumps(file_names))

        else:
            print("Invalid function name")
    else:
        print("Please provide function name and arguments as JSON string")
        ppg_example(savefig=True)