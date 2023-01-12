# An-homogeneous-and-live-database-of-solar-active-regions
This database is a live solar active regions database for solar cycles 23, 24 and 25 and it provides 

# Overview
This python code is the automatic detection method of solar active regions (ARs) of the homogeneous database. The method is based on morphological operation and region growing. It uses synoptic magnetograms from SOHO/MDI and SDO/HMI and they can be got from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/) freely. The detected ARs are given in **allAR.xlsx** file. The time range of ARs is from Carrington rotation (CR) 1909 to CR 2265 now.

Full details of the purpose and design of the code are given in the paper ****.

# Dependencies
The code is tested with Python 3.8.10 on Spyder 5. The only nonstandard library required is astropy. The image processing in our code is based on opencv 4.5.2.

# Usage
First of all, you should download MDI and HMI synoptic magnetograms (**mdi.synoptic_Mr_96m and hmi.Synoptic_Mr_720s**) from [Joint Science Operations Center (JSOC)](http://jsoc.stanford.edu/). The radial magnetograms are used in our code but line-of-sight (LOS) magnetograms can also be processed.

The AR detection modules are in script **ARdetection.py**. Using these modules, there are two ways to process the magnetograms, single processing and batch processing. The script **Single_processing.py** allow you to process a magnetogram once and get the detected ARs and image after each detection module. The **Batch_processing.py** allow you to process all magnetograms once and get an array containing properties of all detected ARs. Before running them, you should first set the path of magnetograms and out path in these two scripts.

The **allAR.xlsx** is the final file of all detected ARs in CR 1909 - CR 2265 (1996-05-05 to 2023-01-01). It provides CR number, label, location of whole AR and each polarity (latitude and longitude), area and flux of each polarity now. More parameters will be given in future.

CR	|Label	|Lat(pos)	|Lon(pos)	|Lat(neg)	|Lon(neg)	|Lat(whole)	|Lon(whole)	|Area(pos, μHem)	|Area(neg, μHem)	|Flux(pos, Mx)|	Flux(neg, Mx)
---|:
1909|	1|	-5.95 |258.77|	-6.93|	247.21|	-6.31	|254.47	1.87E+03	1.63E+03	1.30E+22	-7.71E+21
1909|	2|	11.78 |339.69|	12.84|	336.64|	12.3	|338.21	2.37E+02	2.15E+02	1.07E+21	-1.02E+21
1910|	1|	3.19  | 21.83|	 3.82|  26.66	| 3.49	|24.13	1.57E+02	1.64E+02	7.30E+20	-6.64E+20
1911|	1|	-11.68|	242	 |  -9.02|	239.99|	-10.15|	240.84	4.53E+02	6.79E+02	2.39E+21	-3.24E+21
1911|	2|	12.06	|351.67|	13.53|	357.62|	13.01	| 355.52	2.44E+02	4.39E+02	1.45E+21	-2.66E+21

![image](https://user-images.githubusercontent.com/110174507/212001212-009552ff-1e3b-4011-b147-97a5a33fc4c6.png)


# Note
There is also a script update_latest.py for extending an existing database to a later end date.
According to our calibration, the AR flux of HMI synoptic magnetograms (i.e. ARs after CR 2097(included)) need to be multiplied with a fator 1.36. 

# New options:

# Author
Ruihui Wang, Jie Jiang, Yukun Luo - Beihang University, China
