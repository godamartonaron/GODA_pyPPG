pyPPG example code
==================

In this tutorial you will learn how to use **pyPPG** to engineer morphological PPG biomarkers and export their values.

**Introduction**
----------------
This tutorial provides step-by-step instructions for installing pyPPG and running the example code.

**Step 1**: Install Python 3.10

Download and install Python 3.10 on your computer or server by visiting the official Python website: `Python 3.10 <https://www.python.org/downloads/release/python-3100/>`__.

**Step 2**: Download the Sample PPG Data

You can use the sample PPG data by downloading it from the following link: `Sample PPG data <https://github.com/godamartonaron/GODA_pyPPG/tree/main/sample_data>`__.

**Step 3**: Create and Activate a Virtual Environment

Create a virtual environment named "ppgenv" specifically for Python 3.10 using the py launcher:

.. code-block:: bash

    py -3.10 -m venv ppgenv

Activate the virtual environment:

.. code-block:: bash

   ppgenv\Scripts\activate

**Step 4**: Install pyPPG

While the virtual environment is active, install pyPPG using pip:

.. code-block:: bash

   ppgenv\Scripts\python.exe -m pip install pyPPG

*WARNING*: Running pip as the 'root' user can result in broken permissions and conflicting behaviour with the system package manager. It is recommended to use a virtual environment instead, or please review carefully the `list of package requirements <https://github.com/godamartonaron/GODA_pyPPG/blob/main/docs/requirements.txt>`__.

**Step 5**: Run the Example Code

Open the Python interpreter:

.. code-block:: bash

   python

Run the example code, load the example files (*.mat*, *.txt*, *.csv*, or *.edf* formats) and check the results:

.. code-block:: python

   from pyPPG.example import ppg_example
   ppg_example()

The resulting figures and outcomes are stored within the *temp_dir* folder, which is automatically generated within the project directory.

**Step 6**: Exit the Python Interpreter and Deactivate the Virtual Environment

To exit the Python interpreter, type:

.. code-block:: python

   exit()

Deactivate the virtual environment:

.. code-block:: bash

   deactivate


You have successfully installed pyPPG, executed the example code, and explored the results. Feel free to customize and use pyPPG for your projects.


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
#. **SQI calculation**: This module calculates the PPG Signal Quality Index based on beat template correlation.
#. **Save data**: This module allows for the saving of extracted Fiducial points, Biomarkers, and Statistics into a .csv file.

.. image:: PPG_MAT_sample.png
   :align: center



