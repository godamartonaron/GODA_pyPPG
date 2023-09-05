Comprehensive PPG Analysis
==========================

In this tutorial we will learn how to extract the biomarkers from PPG pulse waves.
Our objectives are to:

    * Detect the standard fiducial points on PPG pulse waves
    * Calculate pulse wave biomarkers from the fiducial points
    * Saving data in different data format

Import Python packages:
________________________

.. code-block:: python

    from pyPPG import PPG, Fiducials, Biomarkers
    from pyPPG.datahandling import load_data, plot_fiducials, save_data
    import pyPPG.fiducials as FP
    import pyPPG.biomarkers as BM
    import pyPPG.ppg_sqi as SQI

    import numpy as np
    import sys
    import json
    import pandas as pd


Setup input parameters:
_______________________

.. code-block:: python

    fs = 100
    start_sig = 0
    end_sig = -1
    correct = True
    correct = True
    savingfolder = 'temp_dir'
    savingformat = 'csv'

Loading a raw PPG signal:
__________________________

.. code-block:: python

    # Load the raw PPG signal
    ppg_data = load_data(fs, start_sig, end_sig, filtering, correct)

    # Create a PPG class
    s = PPG(ppg_data)

Identify fiducial points:
_________________________

.. code-block:: python

    # Init the fiducials package
    fpex = FP.FpCollection(s)

    # Extract fiducial points
    fiducials = fpex.get_fiducials(s, correct=True) + s.start_sig


Plot fiducial points:
_____________________

.. code-block:: python

    # Create a fiducials class
    fp = Fiducials(fiducials)

    # Plot fiducial points
    plot_fiducials(s, fp, savingfolder)

PPG fiducial points
     .. image:: PPG_MAT_sample.png
       :align: center

Calculate PPG biomarkers:
_________________________

.. code-block:: python

    # Init the biomarkers package
    bmex = BM.BmCollection(s, fp)

    # Extract biomarkers
    bm_defs, bm_vals, bm_stats = bmex.get_biomarkers()

    # Create a biomarkers class
    bm = Biomarkers(bm_defs, bm_vals, bm_stats)


Calculate PPG SQI:
_________________________

.. code-block:: python

    # Get PPG SQI
    ppgSQI = round(np.mean(SQI.get_ppgSQI(s.filt_sig, s.fs, fp.sp)) * 100, 2)
    print('Mean PPG SQI: ', ppgSQI, '%')

Save PPG data:
______________

.. code-block:: python

    # Save PPG struct, fiducial points, biomarkers
    save_data(s, fp, bm, savingformat, savingfolder)



Extracted fiducial points
 .. image:: FID_vals.png
   :align: center

Extracted biomarkers
 .. image:: BM_vals.png
   :align: center

Biomarkers statistics
 .. image:: BM_stats.png
   :align: center

Biomarkers definitions
 .. image:: BM_defs.png
   :align: center
