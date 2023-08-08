from DataHandling import*
import FiducialPoints as Fp
import Biomarkers as Bm
from Statistics import*

import matplotlib.pyplot as plt
import numpy as np
from dotmap import DotMap
import time

import sys
import json

###########################################################################
################################## EXAMPLE ################################
###########################################################################
def ppg_example(data_path="",filtering=True,correct=True, savefig=True, savedata=True, savingformat="mat",savingfolder="temp_dir"):
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
    :param filtering: a bool for filtering
    :type filtering: bool
    :param savefig: a bool for current figure saving
    :type savefig: bool
    :param correct: a bool for fiducials points corretion
    :type correct: bool
    :param savedata: a bool for saving fiducial points, biomarkers, and statistics
    :type savedata: bool
    :param savingformat: file format of the saved date, the provided file formats .mat and .csv
    :type savingformat: str
    :param savingfolder: location of the saved data
    :type savingfolder: str

    :return: fiducial points, a dictionary where the key is the name of the fiducial pints and the value is the list of fiducial points
    '''


    ## Loading a raw PPG signal
    s=load_data(data_path,filtering)

    ## Get Fiducial points
    fp = Fp.FiducialPoints(s)
    fiducials=fp.getFiducialPoints(correct)

    if savefig:
        ## Plot Fiducials Points
        plot_fiducials(s, fiducials)

    if savedata:
        ## Get Biomarkers
        bm = Bm.Biomarkers(s, fiducials)
        biomarkers_vals, biomarkers_defs = bm.getBiomarkers()

        ## Get Statistics
        statistics = Statistics(fiducials['sp'], fiducials['on'], biomarkers_vals)

        ## Save data
        save_data(s,fiducials,biomarkers_vals,biomarkers_defs,statistics,savingformat,savingfolder)

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
