import pyPPG
from pyPPG import PPG, Fiducials, Biomarkers
from pyPPG.datahandling import load_data, plot_fiducials, save_data
import pyPPG.preproc as PP
import pyPPG.fiducials as FP
import pyPPG.biomarkers as BM
import pyPPG.ppg_sqi as SQI

import pyPPG.validation.pw_anal as PW

import numpy as np
import sys
import json
import pandas as pd

###########################################################################
################################## EXAMPLE ################################
###########################################################################
def ppg_example(data_path="", fs=0, start_sig=0, end_sig=-1, fiducials=pd.DataFrame(), process_type="both", channel="Pleth",
                filtering=True, fL=0.5000001, fH=12, order=4, sm_wins={'ppg':50,'vpg':10,'apg':10,'jpg':10}, correction=pd.DataFrame(),
                savingfolder="temp_dir", savefig=True, show_fig=True, savingformat="csv", print_flag=True, use_tk=False, check_ppg_len=True):
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
        - ppg: windows for PPG signal
        - vpg: windows for PPG' signal
        - apg: windows for PPG" signal
        - jpg: windows for PPG'" signal
    :type sm_wins: dict
    :param correction: DataFrame where the key is the name of the fiducial points and the value is bool
    :type correction: DataFrame
    :param savingfolder: location of the saved data
    :type savingfolder: str
    :param savefig: a bool for current figure saving
    :type savefig: bool
    :param show_fig: a bool for show figure
    :type show_fig: bool
    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str
    :param print_flag: a bool for print message
    :type print_flag: bool
    :param use_tk: a bool for using tkinter interface
    :type use_tk: bool
    :param check_ppg: a bool for checking ppg length and sampling frequency
    :type check_ppg: bool

    :return:
        - fiducial points: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points
        - s: object of PPG signal

    Example:

        .. code-block:: python

            from pyPPG.example import ppg_example

            # run example code
            ppg_example(savedata=True, savefig=True)

    '''

    ## Loading a raw PPG signal
    signal = load_data(data_path=data_path, fs=fs, start_sig=start_sig, end_sig=end_sig, channel=channel, use_tk=True)

    ## Preprocessing
    # Initialise the filters
    prep = PP.Preprocess(fL=fL, fH=fH, order=order, sm_wins=sm_wins)

    # Filter and calculate the PPG, PPG', PPG", and PPG'" signals
    signal.filtering = filtering
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

    ## PPG SQI
        # Create a fiducials class
        fp = Fiducials(fp=fiducials)

        # Calculate SQI
        ppgSQI = round(np.mean(SQI.get_ppgSQI(ppg=s.ppg, fs=s.fs, annotation=fp.sp)) * 100, 2)
        if print_flag: print('Mean PPG SQI: ', ppgSQI, '%')

    ## Plot fiducial points
        plot_fiducials(s=s, fp=fp, savefig=savefig, savingfolder=savingfolder, show_fig=show_fig, print_flag=print_flag, use_tk=use_tk)

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

    ## Save data
        fp_new = Fiducials(fp=fp.get_fp() + s.start_sig)
        save_data(s=s, fp=fp_new, bm=bm, savingformat=savingformat, savingfolder=savingfolder, print_flag=print_flag)

    if print_flag: print('Program finished')

    return fiducials + s.start_sig, s


###########################################################################
############################## RUN EXAMPLE CODE ###########################
###########################################################################
if __name__ == "__main__":

    if len(sys.argv) > 1:
        input_data = json.loads(sys.argv[1])
        function_name = input_data['function']
        function_args = input_data['args']

        if function_name == 'ppg_example':
            if input_data['args']['process_type'] == 'fiducials':
                fiducials = ppg_example(**function_args)
            if input_data['args']['process_type'] == 'biomarkers':
                fiducials = input_data['args']['fiducials']
                fiducials = fiducials.replace("'", '"')
                fiducials = json.loads(fiducials)
                fiducials = pd.DataFrame(fiducials['data'], columns=fiducials['columns'], index=fiducials['index'])

                def subtract_if_numeric(x):
                    if isinstance(x, (int, float)):
                        tmp_diff = x - input_data['args']['start_sig']
                        return tmp_diff
                    else:
                        return x

                input_data['args']['fiducials'] = fiducials.applymap(subtract_if_numeric)
                function_args = input_data['args']
                ppg_example(**function_args)

            json_data = fiducials.to_json(orient="split")
            print(json.dumps(json_data))

        else:
            print("Invalid function name")
    else:
        print("Please provide function name and arguments as JSON string")
        ppg_example(savefig=True)