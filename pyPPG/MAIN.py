from DataHandling import*
from Prefiltering import*
from FiducialPoints import*
from Biomarkers import*
from Summary import*
from Statistics import*

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
from dotmap import DotMap
from tkinter import filedialog
import mne
import time

from six.moves import cPickle as pickle
import matplotlib.mlab


###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':

    ## Load data
    s=load_data(filtering=1)

    ## Get Fiducials Points
    fiducials = getFiducialsPoints(s,correct=1)

    ## Plot Fiducials Points
    plot_fiducials(s, fiducials)

    ## Get Fiducials Biomarkers, Summary and Statistics
    ppg_biomarkers = Biomarkers(s, fiducials)
    ppg_summary = Summary(s.v, fiducials['sp'], fiducials['on'], s.fs)
    ppg_statistics = Statistics(fiducials['sp'], fiducials['on'], ppg_biomarkers)

    ## Save data
    save_data(fiducials,ppg_biomarkers,ppg_summary,ppg_statistics)

    print('Program finished')
