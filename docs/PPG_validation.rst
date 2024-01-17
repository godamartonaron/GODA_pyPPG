Validation and Benchmarking
============================

.. raw:: html

   <a href="https://colab.research.google.com/drive/1TaaNPaSPj5tq82awgHkBtxV9hXaiE1PT#scrollTo=UXDo7gvgGJEv&uniqifier=6">Colab Notebook</a>


In this tutorial we will learn how to validate and benchmark the fiducial point extraction of photoplethysmogram (PPG) signal.

Our objectives are to:

    * Compare the fiducial points annotations of two manual annotations
    * Validate the fiducial points detection of *pyPPG*
    * Benchmark the fiducial points detection of *pyPPG*, *PulseAnalyse*, *PPGFeat*
    * Make Bland-Altman analysis for the annotations and detection results
    * Saving all results

Download and extract the manual annotations of fiducial points by accessing the provided link: `PPG-BP manual annotations <https://github.com/godamartonaron/GODA_pyPPG/raw/main/pyPPG/validation/PPG-BP_annot.zip>`__.

Additionally, acquire and extract the MATLAB codes essential for evaluating *PPGFeat* and *PulseAnalyse* by using the following link: `MATLAB codes <https://github.com/godamartonaron/GODA_pyPPG/raw/main/pyPPG/validation/MATLAB_codes.zip>`__.

Ensure that you place the extracted contents in the designated folder associated with this script.


Import Python packages:
-----------------------

* Install the pyPPG toolbox for PPG analysis

.. code-block:: python

    pip install pyPPG==1.0.67

* Import required components from pyPPG

.. code-block:: python

    import pyPPG.validation.pw_anal as PW

* Import other packages

.. code-block:: python

    import pandas as pd
    from datetime import datetime
    import os

Initialisation:
-----------------

* Initialise the following input parameters to the validation and benchmarking

.. code-block:: python

    # Initialise the pulse wave package
    pwex = PW.PulseWaveAnal()

    # Initialise the correction
    correction = pd.DataFrame()
    corr_on = ['on', 'v', 'w', 'f']
    corr_off = ['dn', 'dp']
    correction.loc[0, corr_on] = True
    correction.loc[0, corr_off] = False

    # Save time for the results of the output folders
    date = datetime.now()
    dname = str(date.year) + '_' + str(date.month) + '_' + str(date.day) + '_' + str(date.hour) + '_' + str(date.minute)

Validation of the fiducial points:
----------------------------------

* Evaluate the manual annotations and *pyPPG* detection of the fiducial points

.. code-block:: python

    # Run PPG-BP Evaluation for the manual annotation and pyPPG
    pwex.eval_PPG_BP(plts=True, correction=correction, dname=dname, prnt=False)

Extract other fiducial points:
--------------------------------

* Extract the detected the fiducial points of *PPGFeat* and *PulseAnalyse*

.. code-block:: python

    # Command to run MATLAB script for PPGFeat
    current_directory = os.getcwd()
    script_folder = current_directory+os.sep+'PPGFeat'
    scipt='get_PPGFeat_fps'
    pwex.run_matlab_script(script_folder,scipt,dname,'')

    # Command to run MATLAB script for PulseAnalyse
    script_folder = current_directory+os.sep+'PulseAnalyse'
    scipt='get_PA_fps'
    pwex.run_matlab_script(script_folder,scipt,dname,'')

Benchmarking:
--------------

* Compare the results of *PPGFeat* and *PulseAnalyse* with *pyPPG*

.. code-block:: python

    # Run Benchmarking
    pwex.benchmark_PPG_BP(detector='PPGFeat', dname=dname, plt=True, prnt=False)
    pwex.benchmark_PPG_BP(detector='PulseAnal', dname=dname, plt=True, prnt=False)

    # Run Bland-Altman analysis
    script_folder = current_directory+os.sep+'BlandAltman'
    scipt='BlandAltman_anal'
    pwex.run_matlab_script(script_folder, scipt, dname, 'MG_PC')
    pwex.run_matlab_script(script_folder, scipt, dname, 'pyPPG')
    pwex.run_matlab_script(script_folder, scipt, dname, 'PPGFeat')
    pwex.run_matlab_script(script_folder, scipt, dname, 'PulseAnal')

The resulting figures and outcomes are stored within the *results* folder, which is automatically generated within the project directory.

