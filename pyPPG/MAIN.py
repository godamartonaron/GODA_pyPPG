from DataHandling import*
from Prefiltering import*
from FiducialPoints import*
from Biomarkers import*
from Summary import*
from Statistics import*

import matplotlib.pyplot as plt
import numpy as np
from dotmap import DotMap
import time


###########################################################################
####################### Data Acquisition from Files #######################
###########################################################################
if __name__ == '__main__':

    ## Load data
    s=load_data(filtering=True)

    ## Get Fiducials Points
    fiducials = getFiducialPoints(s,correct=True)

    ## Plot Fiducials Points
    plot_fiducials(s, fiducials,savefig=True)

    ## Get Fiducials Biomarkers, Summary and Statistics
    ppg_biomarkers = Biomarkers(s, fiducials)
    ppg_statistics = Statistics(fiducials['sp'], fiducials['on'], ppg_biomarkers)

    ## Save data
    save_data(fiducials,ppg_biomarkers,ppg_statistics)

    print('Program finished')
