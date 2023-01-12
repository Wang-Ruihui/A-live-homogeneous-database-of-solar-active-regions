# An-homogeneous-and-live-database-of-solar-active-regions
This database is a homogeneous and live solar active regions database for solar cycles 23, 24 and 25 and it provides basic parameters of AR. Some useful parameters describing the contribution of ARs to the polar field will be given in future. It can be used to not only the solar cycle variation research but also solar cycle prediction.

# Overview
This python code is the automatic detection method of solar active regions (ARs) of the homogeneous database. The method is based on morphological operation and region growing. It uses synoptic magnetograms from SOHO/MDI and SDO/HMI and they can be got from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/) freely. The detected ARs are given in **allAR.xlsx** file. The time range of ARs is from Carrington rotation (CR) 1909 to CR 2265 at present. ARs in new CR will be included continually.

Full details of the purpose and design of the code are given in the paper ****.

# Dependencies
The code is tested with Python 3.8.10 on Spyder 5. The only nonstandard library required is astropy. The image processing in our code is based on opencv 4.5.2.

# Usage
First of all, you should download MDI and HMI synoptic magnetograms (**mdi.synoptic_Mr_96m and hmi.Synoptic_Mr_720s**) from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/). The radial magnetograms are used in our code but line-of-sight (LOS) magnetograms can also be processed.

The AR detection modules are in script **ARdetection.py**. Using these modules, there are two ways to process the magnetograms, single processing and batch processing. The script **Single_processing.py** allow you to process a magnetogram once and get the detected ARs and image after each detection module. The **Batch_processing.py** allow you to process all magnetograms once and get an array containing properties of all detected ARs. Before running them, you should first set the path of magnetograms and out path in these two scripts.

The **allAR.xlsx** is the final file of all detected ARs in CR 1909 - CR 2265 (1996-05-05 to 2023-01-01). It provides CR number, label, location of whole AR and each polarity (latitude and longitude), area and flux of each polarity at present. More parameters will be given in future.

![image](https://user-images.githubusercontent.com/110174507/212001212-009552ff-1e3b-4011-b147-97a5a33fc4c6.png)


# Note
According to our calibration, the flux of ARs from HMI synoptic magnetograms (i.e. ARs after CR 2097(included)) need to be multiplied with a fator 1.36. 

# Author
Ruihui Wang, Jie Jiang, Yukun Luo - Beihang University, China
