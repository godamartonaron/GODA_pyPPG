pyPPG example code
==================

In this tutorial you will learn how to use **pyPPG** to engineer morphological PPG biomarkers and export their values.

**Introduction**
----------------
Once you have installed the most recent version of Python on your computer, you can proceed to install the pyPPG toolbox by executing the command **pip install pyPPG** in your command line interface. The installation of pyPPG grants access to all the packages within the toolbox. If you encounter any errors during the installation process, please refer to the 'requirements.txt' list and consider upgrading your Python Interpreter accordingly.

**Example pyPPG code**
------------------------
The provided example code consists of seven modules that effectively showcase the capabilities of the pyPPG toolbox.

#. **Raw PPG Signal Loading**: This module facilitates the loading of raw PPG signals from various file formats, including .mat, .csv, .txt, or .edf.
#. **Fiducial Point Extraction**: This module focuses on extracting fiducial points from PPG signals, encompassing PPG, PPG', PPG'', and PPG'''.
#. **Fiducial Points Plotting**: Here, the extracted fiducial points are visually represented through plotting.
#. **Biomarker Extraction**: This module offers the extraction of 74 distinct PPG biomarkers, categorized into:

    I. PPG signal characteristics
    II. Signal ratios
    III. PPG derivatives
    IV. Derivative ratios

#. **Biomarker Statistics**: A concise summary of the 74 PPG biomarkers is provided within this module.
#. **SQI calculation**: This module calculates PPG Signal Quality Index based on beat template correlation.
#. **Save data**: This module allows for the saving of extracted Fiducial points, Biomarkers, and Statistics into a .csv file.

The resulting figures and outcomes are stored within the *temp_dir* folder, which is automatically generated within the project directory.

.. image:: PPG_MAT_sample.png
   :align: center



