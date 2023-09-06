from pyPPG import PPG, Fiducials, Biomarkers
from pyPPG.preproc import Preprocessing
from pyPPG.datahandling import load_data, plot_fiducials, save_data
import pyPPG.fiducials as FP
import pyPPG.biomarkers as BM
import pyPPG.ppg_sqi as SQI

import numpy as np
import sys
import json
import pandas as pd


###########################################################################
################################## EXAMPLE ################################
###########################################################################
def ppg_example(data_path="", fs=[], start_sig=0, end_sig=-1, filtering=True, correct=True, process_type="both",
                savingfolder="temp_dir", savefig=True, savingformat="csv", fiducials=[], print_flag = False):
    '''
    This is an example code for PPG analysis. The main parts:
        1) Loading a raw PPG signal: various file formats such as .mat, .csv, .txt, or .edf.
        2) Get Fiducial points: extract the fiducial points of PPG, PPG', PPG'' and PPG'" signals
        3) Plot Fiducials Points
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
    :param start_sig: beginning the of signal in sample
    :type start_sig: int
    :param fs: sampling_frequency
    :type fs: int
    :param end_sig: end of the signal in sample
    :type end_sig: int
    :param filtering: a bool for filtering
    :type filtering: bool
    :param correct: a bool for fiducials points corretion
    :type correct: bool
    :param process_type: the type of the process, which can be "fiducials", "biomarkers", or "both"
    :type process_type: str
    :param savingfolder: location of the saved data
    :type savingfolder: str
    :param savefig: a bool for current figure saving
    :type savefig: bool
    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str
    :param print_flag: a bool for print message
    :type print_flag: bool

    :return: fiducial points: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points

    Example:

        .. code-block:: python

            from pyPPG.example import ppg_example

            # run example code
            ppg_example(savedata=True, savefig=True)

    '''

    ## Loading a raw PPG signal
    ppg_data = load_data(data_path, fs, start_sig, end_sig)

    ## Preprocessing
    if filtering:
        ppg_data.filt_sig, ppg_data.filt_d1, ppg_data.filt_d2, ppg_data.filt_d3 = Preprocessing(ppg_data, filtering=filtering)

    ## Create a PPG class
    ppg_data.filtering = filtering
    ppg_data.correct = correct
    s = PPG(ppg_data)

    ## Get Fiducial points
    if process_type == 'fiducials' or process_type == 'both':
        # Init the fiducials package
        fpex = FP.FpCollection(s)

        # Extract fiducial points
        fiducials = fpex.get_fiducials(s, correct) + s.start_sig
        if print_flag: print("Fiducial points:\n", fiducials)

    if savefig:
        # Create a fiducials class
        fp = Fiducials(fiducials)

        # Plot fiducial points
        plot_fiducials(s, fp, savingfolder, print_flag)

    ## Get Biomarkers and Statistics
    if process_type == 'biomarkers' or process_type == 'both':
        # Init the biomarkers package
        fp = Fiducials(fiducials)
        bmex = BM.BmCollection(s, fp)

        # Extract biomarkers
        bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()

        if print_flag:
            tmp_keys = bm_stats.keys()
            print('Statistics of the biomarkers:')
            for i in tmp_keys: print(i, '\n', bm_stats[i])

        # Create a biomarkers class
        bm = Biomarkers(bm_defs, bm_vals, bm_stats)

        ## Save data
        save_data(s, fp, bm, savingformat, savingfolder, print_flag)

    ## PPG SQI
    fp = Fiducials(fiducials)
    ppgSQI = round(np.mean(SQI.get_ppgSQI(s.filt_sig, s.fs, fp.sp)) * 100, 2)
    if print_flag: print('Mean PPG SQI: ', ppgSQI, '%')

    if print_flag: print('Program finished')
    return fiducials


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