from pyPPG import*
from datahandling import*
import pyPPG.fiducials as FP
import pyPPG.biomarkers as BM

import sys
import json

###########################################################################
################################## EXAMPLE ################################
###########################################################################
def ppg_example(data_path="",start = 0, end = 0, filtering=True, correct=True, savingfolder="temp_dir", savefig=True, savedata=True, savingformat="csv"):
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
        6) Save data: save the extracted Fiducial points, Biomarkers, and Statistics into .csv file

    :param data_path: path of the PPG signal
    :type data_path: str
    :param start: beginning the of signal in sample
    :type start: int
    :param end: end of the signal in sample
    :type end: int
    :param filtering: a bool for filtering
    :type filtering: bool
    :param correct: a bool for fiducials points corretion
    :type correct: bool
    :param savingfolder: location of the saved data
    :type savingfolder: str
    :param savefig: a bool for current figure saving
    :type savefig: bool
    :param savedata: a bool for saving fiducial points, biomarkers, and statistics
    :type savedata: bool
    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str

    :return: fiducial points, a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points


    Example:

        .. code-block:: python

            from pyPPG import ppg_example

            # run example code
            ppg_example(savedata=True, savefig=True)

    '''

    ## Loading a raw PPG signal
    ppg_data = load_data(data_path,start,end,filtering)
    s = PPG(ppg_data)

    ## Get Fiducial points
    fpex = FP.FpCollection(s)
    fiducials=fpex.get_fiducials(s,correct)+s.start
    fp = Fiducials(fiducials)

    if savefig:
        ## Plot Fiducials Points
        plot_fiducials(s, fp, savingfolder)

    if savedata:
        ## Get Biomarkers and Statistics
        bmex = BM.BmCollection(s, fp)
        bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()
        bm = Biomarkers(bm_defs, bm_vals , bm_stats)

        ## Save data
        save_data(s, fp, bm, savingformat, savingfolder)

    print('Program finished')
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
            fiducials = ppg_example(**function_args)
            json_data = fiducials.to_json(orient="split")
            print(json.dumps(json_data))
        else:
            print("Invalid function name")
    else:
        print("Please provide function name and arguments as JSON string")
        ppg_example(savedata=True, savefig=True)