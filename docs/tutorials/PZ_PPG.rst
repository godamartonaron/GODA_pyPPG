PhysioZoo PPG morphological analysis
====================================

In this tutorial you will learn how to use **PhysioZoo PPG** to calculate morphological PPG biomarkers (i.e. pulse wave features) and export their values.

**Introduction**
----------------------
The PPG signal is an optical measurement of the arterial pulse wave `(Charlton et
al. 2019) <https://journals.physiology.org/doi/full/10.1152/ajpheart.00218.2019>`__, *i.e.*, the wave generated when blood is ejected from the heart, temporarily increasing arterial pressure and causing vessel expansion and contraction `(Alastruey et
al. 2023) <https://journals.physiology.org/doi/full/10.1152/ajpheart.00705.2022>`__, the PPG signal is influenced by a range of physiological systems, such as: the heart, including heart rate, heart rhythm, and the nature of ejection `(Charlton et
al. 2020) <https://ieeexplore.ieee.org/abstract/document/9733047/>`__; the blood vessels, including vessel stiffness, diameter, and blood pressure; the microvasculature, including peripheral compliance and resistance `(Charlton et
al. 2020) <https://ieeexplore.ieee.org/abstract/document/9733047/>`__; the autonomic nervous system which influences heart rate variability `(Gil et
al. 2010) <https://iopscience.iop.org/article/10.1088/0967-3334/31/9/015/meta>`__; and the respiratory system, which impacts the pulse wave through changes in intrathoracic pressure `(Charlton et
al. 2017) <https://iopscience.iop.org/article/10.1088/1361-6579/aa670e/meta>`__. Thus, there is potential to extract much physiological information from the PPG signal.

Studying the morphological characteristics of the PPG may provide information on cardiovascular health.
**PhysioZoo PPG** provides a framework and tools for extracting morphological biomarkers from the PPG signal.

**Performing PPG morphological analysis**
------------------------------------------------------------
Start by entering the PPG interface by clicking on the 'Pulse' menu on the top left, then load some PPG example by clicking File -> Open data file -> ppg_example.txt. The program will automaticly present the PPG file you imported.

.. .. image:: before_analysis.png
   :align: center

To perform the analysis, please follow the instructions:

#. Prefiltering the signal: On the left panel, select the "Configuration" tab. On the bottom of the tab, you will find a section labeled: **Fiducials filtering parameters**. The following
   filters have been implemented as default in the pyPPG toolbox:

    * **Bandpass filtering between 0.5-12 Hz**: A fourth-order Chebyshev Type II filter was used for the original signal. The 12 Hz low-pass cut-off was used to avoid time-shifting of fiducial points (particularly pulse onset, and dicrotic notch) and to eliminate unwanted high-frequency content from the PPG derivatives. The 0.5 Hz high-pass cut-off was used to minimize baseline wandering whilst retaining content at low heart rates.
    * **20 ms moving average filtering (MAF)**: In the case of very noisy signals, some high-frequency content can remain in the band-pass filter signal. For this purpose, a 20 ms standard flat (boxcar or top-hat) MAF with a 22.5 Hz cut-off frequency was applied after the band-pass filtering.
    * **10 ms MAF for the PPG derivatives**: To eliminate the high-frequency content in the PPG derivatives, a 10 ms standard flat (boxcar or top-hat) MAF with 45 Hz cut-off frequency was applied.

#. Definition of the window for anlysis: On the right panel, define the W.S. (start of the window) and the W.L. (length of the window) you want to analyze. You can analyze all of your signal or part of it. Note that if you analyze a long window,it may take some time.

#. Click the **Find Fiducials** button. The fiducial poits will be detected and highlighted while the biomarkers will be automatically engineered and displayed on the lower pannels.

Congrats! You have made your first morphological analysis with **PhysioZoo PPG**!
The biomarkers are divided into two different categories: Duration and Amplitudes, the statistical measurments of the biomarkers will be presented in a table, in the bottom panel.

.. .. image:: after_analysis.png
   :align: center

.. note:: For PPG anlysis 9 statistical mesurment computed over the selected window (defined by W.S. W.L.) will be presented for each biomarker namely: signal duration; average (AVG); median (MED); standard deviation (STD); lower and upper quartiles (Q1, Q3); inter-quartile range (IQR); Skewness (SKW, indicating a lack of symmetry in the distribution; Kurtosis (KUR, indicating the pointedness of a peak in the distribution curve); and the average difference between the mean and each data value (MAD)

**Exporting fiducial points**
--------------------------------------------

You can export the fiducial points. Go to File -> Save fiducial points. The excel file contains the computed fiducial points for each lead.


.. .. image:: results_fiducials.png
   :align: center

**Exporting morphological biomarkers**
--------------------------------------------

You can export the morphological biomarkers. Go to File -> Save fiducial biomarkers. The excel file contains the engineered PPG biomarkers.

.. .. image:: results_mor_analysis.PNG
   :align: center