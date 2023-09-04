from pyPPG import PPG, Fiducials, Biomarkers
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
def ppg_example(data_path="", start_sig=0, fs=[], end_sig=0, filtering=True, correct=True,
                process_type="both",savingfolder="temp_dir", savefig=True, savingformat="csv", fiducials=[]):
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

    :return: fiducial points: DataFrame where the key is the name of the fiducial pints and the value is the list of fiducial points

    Example:

        .. code-block:: python

            from pyPPG.example import ppg_example

            # run example code
            ppg_example(savedata=True, savefig=True)

    '''

    ## Loading a raw PPG signal
    ppg_data = load_data(data_path, fs, start_sig, end_sig, filtering, correct)
    s = PPG(ppg_data)

    if process_type == 'fiducials' or process_type == 'both':
        ## Get Fiducial points
        fpex = FP.FpCollection(s)
        fiducials = fpex.get_fiducials(s, correct) + s.start_sig

    if savefig:
        ## Plot Fiducials Points
        fp = Fiducials(fiducials)
        plot_fiducials(s, fp, savingfolder)

    if process_type == 'biomarkers' or process_type == 'both':
        ## Get Biomarkers and Statistics
        fp = Fiducials(fiducials)
        bmex = BM.BmCollection(s, fp)
        bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()
        bm = Biomarkers(bm_defs, bm_vals, bm_stats)

        ## Save data
        save_data(s, fp, bm, savingformat, savingfolder)

    # PPG SQI
    fp = Fiducials(fiducials)
    ppgSQI = round(np.mean(SQI.get_ppgSQI(s.filt_sig, s.fs, fp.sp)) * 100, 2)
    print('Mean PPG SQI: ', ppgSQI, '%')

    # print('Program finished')
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