# Active Region database for Influence on Solar cycle Evolution (ARISE)
# (formerly a live homogeneous database of solar active regions)

# Overview
This is the `**A**ctive **R**egion database for **I**nfluence on **S**olar cycle **E**volution (ARISE)', developed by **Ruihui Wang, Jie Jiang, and Yukun Luo**, Beihang University, China. 

We offer two versions of the database: one includes all detected ARs (**database_allAR.xlsx**), and the other excludes repeat ARs (**database_remove_repeat_AR.xlsx**). Both versions contain the same parameters, including basic parameters and those related to solar cycle variability (dipole fields, calculated with $\lambda_R = 6.44^\circ$). However, caution should be taken when using the dipole field data from _database_allAR.xlsx_, as some long-lived ARs are recorded multiple times.

The database covers ARs from Carrington Rotations (CRs) 1909 to 2303, corresponding to the period from May 1996 to October 2025, encompassing cycles 23, 24, and 25. Data from subsequent CRs will be continuously added.

For a detailed explanation of the purpose and design of the database, please refer to the following publications:

1. Ruihui Wang, Jie Jiang, Yukun Luo, "Toward a Live Homogeneous Database of Solar Active Regions Based on SOHO/MDI and SDO/HMI Synoptic Magnetograms. I. Automatic Detection and Calibration," ApJS, https://doi.org/10.3847/1538-4365/acef1b.
2. Ruihui Wang, Jie Jiang, Yukun Luo, "II. Parameters for Solar Cycle Variability," ApJ, https://doi.org/10.3847/1538-4357/ad5b5f.

# Dependencies
The code has been tested with Python 3.8.10 on Spyder 5. The nonstandard libraries required are _astropy_ and _sunpy._ Image processing in the code relies on _opencv_ version 4.5.5

# Usage
These Python codes include the automatic detection method for solar active regions (ARs) and the programs for creating the database. The detection method is based on morphological operations and region growing. We use synoptic magnetograms from SOHO/MDI and SDO/HMI, which can be freely obtained from the [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/).

To use the codes, you first need to download MDI and HMI synoptic magnetograms (_mdi.synoptic_Mr_96m and hmi.Synoptic_Mr_720s_) from JSOC. The radial magnetograms are used in our code, but line-of-sight (LOS) magnetograms can also be processed with appropriate adjustments to the detection parameters


**ARdetection.py** contains the AR detection methods, including the original detection method (**ARdetection_org**) and the updated method with repeat-AR-removal (**ARdetection**). Our detection process consists of five modules: adaptive intensity threshold segmentation, morphological closing and opening operations, region growing, small region removal, and unipolar region removal. The repeat-AR-removal module is integrated between the small region removal and unipolar region removal steps in the **ARdetection** function.

Additionally, we provide two functions in the file to process a map: **Get_ARi**, which generates the detection image, and **Get_ARP**, which returns the parameters of the detected ARs.

**ARparameters.py** describes the methods to calculate the parameters of each detected AR, such as area, flux, and dipole fields.

**Single_processing_MG.py** enables you to process a single magnetogram and obtain the detected ARs labeled on the map.

**Batch_processing_MG_param.py** enables you to process all magnetograms at once and generates an array containing the parameters of all detected ARs.

**database_allAR.xlsx** is the file of all detected ARs in CR 1909 - CR 22303 (1996-05 to 2025-10). It provides the CR number, label, the latitude and longitude of the flux-weighted centroid for both polarities and the entire AR, the area, the flux of each polarity, the maximum magnetic field of the entire AR, the dipole fields, and dipole fields with BMR approximations. Dipole fields are calculated with transport parameter $\lambda_R = 6.44^\circ$. More parameters will be given in the future if necessary. The following picture is a part of the file.
<img width="2987" height="973" alt="image" src="https://github.com/user-attachments/assets/236800c7-cfd0-45c6-96e5-607c0f16e1b7" />

**database_remove_repeat_AR.xlsx** is the same as **database_allAR.xlsx**, except that repeat ARs have been removed from this file.

**Output_ARmap.py** is used to generate maps of each detected AR. The AR maps are named according to their emergence time relative to the specified time0 (default: CR 1912, 1996-07-25 21:35:44.193, which can be changed in codes). Each map is a single AR. ARs emerging on the same day are labeled as *day_00115_1* and *day_00115_2*.

**Output_BoundaryARs.py** is used to detect ARs on the boundary of synoptic maps and generate their maps. Others are the same as **Output_ARmap.py**.

**Batch_processing_MG_outputAR.py** is used to output maps of the detected ARs across all synoptic magnetograms.

**Smooth_Resize_img.py** is a function file that contains functions to change the coordinate of the map, smooth the map, and resize the map to a lower resolution.

**ARmaps.zip** is the outputted low-resolution (180*360) maps of detected AR. The AR maps are named according to their emergence time relative to the specified time0 (default: CR 1912, 1996-07-25 21:35:44.193). 

**ARmaps(Luo2026).zip** provides the AR maps used in Luo, Jiang, and Wang (2026). The detected ARs from the MDI synoptic maps show slight differences from those in the **ARmaps.zip**. 

# Note
According to our calibration, the fluxes, Bmax, and dipole field of ARs from HMI synoptic magnetograms (i.e., ARs after CR 2097(included)) need to be multiplied by a factor of 1.36. 

# Acknowledgements
The database is supported by SCOSTEP/PRESTO 2024 database grant, the National Natural Science Foundation of China Nos. 12173005 and 12350004, and the National Key R&D Program of China No. 2022YFF0503800.

