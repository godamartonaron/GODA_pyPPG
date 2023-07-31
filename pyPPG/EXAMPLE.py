from DataHandling import*
from Preprocessing import*
import FiducialPoints as Fp
import Biomarkers as Bm
from Statistics import*

import matplotlib.pyplot as plt
import numpy as np
from dotmap import DotMap
import time

import sys

###########################################################################
################################## EXAMPLE ################################
###########################################################################
def example_code(filtering=True,correct=True,savefig=True):
    '''
    This is an example code for PPG analysis. The main parts:
        1) Loading a raw PPG signal: various file formats such as .mat, .csv, .txt, or .edf.
        2) Get Fiducial points: extract the fiducial points of PPG, PPG', PPG'' and PPG'" signals
        3) Plot Fiducials Points
        4) Get Biomarkers: extract 74 PPG biomarkers in four categories
            - PPG signal
            - Signal ratios
            - PPG derivatives
            - Derivatives ratios
        5) Get Statistics: summary of the 74 PPG biomarkers
        6) Save data: save the extracted Fiducial points, Biomarkers, and Statistics into .csv file

    :param filtering: a bool for filtering
    :param savefig: a bool for fiducial points saving
    :param correct: a bool for fiducials points corretion
    '''

    python_executable_path = sys.executable

    ## Loading a raw PPG signal
    s=load_data(filtering)

    ## Get Fiducial points
    fp = Fp.FiducialPoints(s)
    fiducials=fp.getFiducialPoints(correct)

    ## Plot Fiducials Points
    plot_fiducials(s, fiducials,savefig)

    ## Get Biomarkers
    bm = Bm.Biomarkers(s, fiducials)
    biomarkers = bm.getBiomarkers()

    ## Get Statistics
    statistics = Statistics(fiducials['sp'], fiducials['on'], biomarkers)

    ## Save data
    save_data(s,fiducials,biomarkers,statistics)

    print('Program finished')



###########################################################################
############################## RUN EXAMPLE CODE ###########################
###########################################################################
example_code()