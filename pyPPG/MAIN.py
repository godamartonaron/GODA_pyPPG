from DataHandling import*
from Prefiltering import*
import FiducialPoints as Fp
import Biomarkers as Bm
from Statistics import*

import matplotlib.pyplot as plt
import numpy as np
from dotmap import DotMap
import time

###########################################################################
################################### MAIN ##################################
###########################################################################
if __name__ == '__main__':

    ## Load data
    s=load_data(filtering=True)

    ## Get Fiducials Points
    fp = Fp.FiducialPoints(s)
    fiducials=fp.getFiducialPoints(correct=True)

    ## Plot Fiducials Points
    plot_fiducials(s, fiducials,savefig=True)

    ## Get Fiducials Biomarkers
    bm = Bm.Biomarkers(s, fiducials)
    biomarkers = bm.getBiomarkers()

    ## Get Statistics
    statistics = Statistics(fiducials['sp'], fiducials['on'], biomarkers)

    ## Save data
    save_data(s,fiducials,biomarkers,statistics)

    print('Program finished')
