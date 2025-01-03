# An Active Region Database for Solar Cycle Variability and Prediction 
# (formerly a live homogeneous database of solar active regions)

# Overview
This is `An Active Region Database for Solar Cycle Variability and Prediction', developed by Ruihui Wang, Jie Jiang, and Yukun Luo. We offer two versions of the database: one includes all detected ARs (**database_allAR.xlsx**), and the other excludes repeat ARs (**database_remove_repeat_AR.xlsx**). Both versions contain the same parameters, including basic parameters and parameters for solar cycle variability (dipole fields). However, caution should be exercised when using dipole field data due to the repeated records of some long-lived ARs.

The temporal coverage of the ARs spans from Carrington Rotations (CRs) 1909-2290, corresponding to the period from May 1996 to November 2024, covering cycles 23, 24, and 25. Data from new CRs will be continuously added.

For a detailed explanation of the purpose and design of the database, please refer to the following publications:

1. "Toward a Live Homogeneous Database of Solar Active Regions Based on SOHO/MDI and SDO/HMI Synoptic Magnetograms. I. Automatic Detection and Calibration," ApJS, https://doi.org/10.3847/1538-4365/acef1b.
2. "II. Parameters for Solar Cycle Variability," ApJ, https://doi.org/10.3847/1538-4357/ad5b5f.
    
We offer two versions of the database: one includes all detected ARs (**database_allAR.xlsx**), and the other excludes repeat ARs (**database_remove_repeat_AR.xlsx**). However, caution should be exercised when using dipole field data due to the repeated records of some long-lived ARs. 

This database is a live homogeneous database of solar active regions for solar cycles 23, 24, and 25.  It can be used for not only the active region long-time variation research (space climate) but also solar cycle prediction.

This python code is the automatic detection method of solar active regions (ARs) of the homogeneous database. The method is based on morphological operation and region growing. It uses synoptic magnetograms from SOHO/MDI and SDO/HMI, which an be obtained from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/) freely. The detected ARs are given in **allAR.xlsx** file. The time range of ARs is from Carrington rotation (CR) 1909 to CR 2271 at present. ARs in the new CR will be included continually.

Full details of the purpose and design of the code are given in the paper [Ruihui Wang, Jie Jiang, Yukun Luo, Toward a live homogeneous database of solar active regions based on SOHO/MDI and SDO/HMI synoptic magnetograms. I. Automatic detection and calibration (APJS)](https://doi.org/10.3847/1538-4365/acef1b).

# Dependencies
The code is tested with Python 3.8.10 on Spyder 5. The nonstandard libraries required are astropy, and sunpy. The image processing in our code is based on opencV 4.5.5.

# Usage
First of all, you should download MDI and HMI synoptic magnetograms (**mdi.synoptic_Mr_96m and hmi.Synoptic_Mr_720s**) from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/). The radial magnetograms are used in our code but line-of-sight (LOS) magnetograms can also be processed with proper changes to the detection parameters.

**ARdetection.py** contains the five AR detection modules, adaptive intensity threshold segmentation, morphological closing operation, and opening operation, region growing, small regions removal, and unipolar regions removal. They can be used by function **ARdetection(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4)**. The meaning and value of function parameters are as follows:

    img_input: the synoptic magnetogram input;    
    # closing operation kernel in module2
    Kernel1 = np.ones((3, 3))
    # opening operation kernel in module2
    Kernel2 = np.ones((11, 11))
    # region growing threshold
    Thresh = 50
    # closing operation kernel in module4
    Kernel3 = np.ones((5, 5))
    # ARs area threshold in module4
    Size = 351
    # dilating operation kernel in module5
    Kernel4 = np.ones((23, 23))
The value of *Kernel2* and *Thresh* is just for MDI magnetograms. For HMI magnetograms, 
    
    Kernel2 = np.ones((9, 9))
    Thresh = 30

Using the function, there are two ways to process the magnetograms, single processing and batch processing. The script **Single_processing.py** allows you to process a magnetogram once and gets the ARs of the map and the image after each detection module. The **Batch_processing.py** allows you to process all magnetograms once and gets an array containing properties of all detected ARs. Before running them, you should first set the path of magnetograms and output path in these two scripts. 

The **ARparameters.py** describes the methods to calculate the location, area, and flux of each detected AR.

The **allAR.xlsx** is the final file of all detected ARs in CR 1909 - CR 2271 (1996-05-05 to 2023-06-14). It provides the CR number, label, latitude and longitude of each polarity and the whole AR, area, and flux of each polarity at present. More parameters will be given in the future. The following picture is a part of the file.

![image](https://user-images.githubusercontent.com/110174507/212001212-009552ff-1e3b-4011-b147-97a5a33fc4c6.png)

# Note
According to our calibration, the fluxes of ARs from HMI synoptic magnetograms (i.e. ARs after CR 2097(included)) need to be multiplied with a factor of 1.36. 

# Author
Ruihui Wang, Jie Jiang, Yukun Luo - Beihang University, China
