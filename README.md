# An-homogeneous-and-live-database-of-solar-active-regions
This database is a live solar active regions database for solar cycles 23, 24 and 25.

# Overview
This python code is the automatic detection method of solar active regions (ARs) of our homogeneous database. The method is based on morphological operation and region growing. It uses synoptic magnetograms from SOHO/MDI and SDO/HMI and they can be got from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/) freely.

This Python code generates a database of bipolar magnetic regions by automated analysis of the Spaceweather HMI Active Region Patch (SHARP) data. The original data are pulled in automatically using sunpy and are provided by the Helioseismic and Magnetic Imager on Solar Dynamics Observatory.

Full details of the purpose and design of the code are given in the paper ****.

# Dependencies
The code is tested with Python 3.8.10 on Spyder 5. The only nonstandard libraries required are astropy.

# Usage
Run the script main.py. You should first set the parameters in this script, especially the time ranges and the output path.

Various output files are produced at each stage of the script:

a file allsharps.txt with a single entry for each detected SHARP, with parameters of the fitted BMR for all SHARPs deemed suitable for BMR fitting (i.e. with good=1 - see the paper for details).
.png images in the directory outputpath showing all of the SHARPs, and the fitted BMRs where appropriate (classified as good or bad depending on whether a BMR is fitted or not).
netcdf files containing the magnetic field of each good SHARP, for use in the next stage below (or to drive other numerical simulations).
a file repeatpairs.txt listing pairs of BMRs identified as repeats.
.png images showing the time evolution of the axial dipole moment predicted by the SFT model for both the original SHARP and the fitted BMR.
a final database file bmrsharps_evol.txt listing all good SHARPs from allsharps.txt with columns for the predicted asymptotic dipole moment added (both for the SHARP and the BMR). Here is an example:
SHARPs from 2010-07-01 00:00:00 to 2010-09-01 00:00:00
-- Produced by anthony.yeates[at]durham.ac.uk --
11
Grid resolution: 180 x 360, smoothing_param = 4, magtype = los, method = cm, maxlon = 90
Selection criteria: (i) sep > 1 deg,  (ii) |imbalance| <  0.5
Last two columns use 10-year 1D SFT simulation with eta=350 km^2/s, v0=0.015 km/s, p=2.33, no decay term.

Be aware that the code will take some time to run if you are trying to cover long time periods. However, if it fails at any point during data download, then you should be able to run it again and it will continue from where it left off. You can force it to start from scratch by setting restart=False.

By default the code is set up to run for a two-month period, and the corresponding allsharps.txt, repeatpairs.txt and bmrsharps_evol.txt files are provided here for verification.

# Note
There is also a script update_latest.py for extending an existing database to a later end date.
According to our calibration, the AR flux of HMI synoptic magnetograms (i.e. ARs after Carrington rotation (CR) 2097(included) ) need to multiply a fator 1.36. 

# New options:

# Author
Ruihui Wang, Jie Jiang, Yukun Luo - Beihang University, China
