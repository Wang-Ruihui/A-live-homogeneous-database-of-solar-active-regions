# -*- coding: utf-8 -*-
"""
process the single magnetogram
"""
import numpy as np
import matplotlib.pyplot as plt
import cv2
import scipy
import pandas as pd
import time
from sunpy.coordinates.sun import carrington_rotation_time
from ARdetection import imshow, Get_ARi
from ARparameters import ARArea, ARFlux, ARLocat, Distance, FDF, IDF, FDFrs, FDFBMR, ARBmax
# ------------------------------------------------------------------------------
# CR1968, CR2155 synoptic maps


def ARparam(img_label, img_input, lamr=6.5):
    # give the AR parameter
    # transport parameters
    num = np.max(img_label)
    CR_array = np.ones(num)*cr
    area = ARArea(img_label, img_input)
    flux = ARFlux(img_label, img_input)
    location = ARLocat(img_label, img_input)
    # the max unsigned magnetic field of AR
    Bmax = ARBmax(img_label, img_input)

    idf = IDF(img_label, img_input)
    fdf = FDF(img_label, img_input, lamr)
    # the BMR approximation of idf and fdf
    idfbmr, fdfbmr = FDFBMR(flux, location, lamr)
    # fdfsm = FDFsm(img_label, img_input, lamr)

    param = np.zeros((17, num))
    param[0, :] = CR_array
    param[1, :] = np.linspace(1, num, num)
    param[2:8, :] = location
    param[8:10, :] = area
    param[10:12, :] = flux
    param[12, :] = Bmax
    param[13, :] = idf
    param[14, :] = fdf
    param[15, :] = idfbmr
    param[16, :] = fdfbmr

    return num, location, param


def match_noaa_to_ar(df2, Locat, cr, img_shape, contour_list, i):
    """

    Parameters
    ----------
    df2 : pd.DataFrame
        dataframe of selected NOAA sunspot, closest to the CM
    Locat : array,list
        list of detected ARs location
    cr : int
        carrington rotation
    img_shape : TYPE
        shape of img
    contour_list : TYPE
        DESCRIPTION.
    i : int
        image label

    Returns
    -------
    None.

    """
    lon = Locat[5, i]

    t1 = carrington_rotation_time(cr)
    t2 = carrington_rotation_time(cr + 1)
    tb = t2 + lon / 360 * (t1 - t2)
    tb = tb.value
    year = int(tb[0:4])
    month = int(tb[5:7])
    day = int(tb[8:10])

    time_mask = (
        (df2['Year'] == year) &
        (df2['Month'] == month) &
        (df2['Day'] >= day - 4) &
        (df2['Day'] <= day + 4)
    )

    if not time_mask.any():
        return None

    matched_rows = df2.loc[time_mask]
    matched_rows = matched_rows[abs(matched_rows['CMD']) < 60]

    if matched_rows.empty:
        return None

    noaa_list = []
    for _, row in matched_rows.iterrows():
        noaa_lon = row['Lon']
        noaa_lat = row['Lat']

        lat_index = (np.sin(noaa_lat / 180 * np.pi) + np.sin(np.pi / 3)) \
            / (2 * np.sin(np.pi / 3)) * img_shape[0]
        lon_index = noaa_lon / 360 * img_shape[1]

        x = float(lon_index)
        y = float(lat_index)

        for contour in contour_list[i]:
            result = cv2.pointPolygonTest(contour, (x, y), measureDist=True)
            # if the distance between a NOAA sunspot and AR smaller than 25 pixel
            if result >= -25:
                noaa_list.append(int(row['NOAA No']))
                break

    if len(noaa_list) == 1:
        return noaa_list[0]
    elif noaa_list:
        return noaa_list
    else:
        return None


t_s = time.time()

cr = 2301

# read fits file and get the synoptic magnetogram data
img_input, module1, module2, module3, module4, module42, module5, \
    result, img_label, img_border = \
    Get_ARi(cr, allregion=False, RemoveRepeat=True)

# results:img_input, module1, module2, module3, module4,module5,  result, img_label, img_border

# get AR parameters
lamr = 6.44  # 7.17 # pre-set transport parameters
num, location, param = ARparam(img_label, img_input, lamr)

# calculate the contour of each AR
ar_contours = []
for i in range(num):
    ar = np.zeros_like(img_label)
    ar[img_label == (i + 1)] = 255
    ar = ar.astype('uint8')
    ar_closed = cv2.morphologyEx(ar, cv2.MORPH_CLOSE, np.ones((101, 101)))
    contours, _ = cv2.findContours(
        ar_closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    ar_contours.append(contours)


# give the NOAA number to each AR
# the results seem good; but need more test
File_CM = 'D:/python program/活动区识别/data/RGO and USAFNOAA/sunspot_CM.xlsx'
df2 = pd.read_excel(File_CM, sheet_name=0)

ar_noaanum = []
for i in range(num):
    noaa_num = match_noaa_to_ar(
        df2, location, cr, img_label.shape, ar_contours, i)
    ar_noaanum.append(noaa_num)
# ------------------------------------------------------------------------------
'''
# figure1: flowchart of the detection modules
plt.rc('font', family='Times New Roman')
plt.rcParams["font.weight"] = 'bold'
plt.rcParams["axes.labelweight"] = 'bold'
plt.rcParams["mathtext.fontset"] = 'stix'

Y, X = img_input.shape
nl = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)','(h)']

fig, axes = plt.subplots(8, 1,  figsize=(6, 9.5))
img = [img_input, module1, module2,
       module3, module4, module42, module5, result]

for i in range(len(axes.flat)):
    ax = axes.flat[i]
    if (i > 0) and i < 7 :
        imshow(ax, img[i],vmax=1,vmin=0)
    else:
        imshow(ax, img[i])
    ax.text(3700, 600, nl[i], fontsize=24)
    if i == 0:
        ax.contour(img_border, colors='darkorange')

axes[0].tick_params(length=4, width=1, labelsize=18)
axes[0].set_xticks([])
axes[0].set_yticks(np.linspace(0, Y, 3))
axes[0].set_yticklabels(
    ['sin(-60$^\circ$)', '0', 'sin(60$^\circ$)'], fontsize=14)
axes[6].tick_params(length=4, width=1, labelsize=18)
axes[6].set_xticks(np.linspace(0, X, 5))
axes[6].set_xticklabels(['0', '90', '180', '270', '360'])
axes[6].set_yticks(np.linspace(0, Y, 3))
axes[6].set_yticklabels(
    ['sin(-60$^\circ$)', '0', 'sin(60$^\circ$)'], fontsize=14)
# remove ticks
for i in range(1, len(axes.flat)-1):
    ax = axes.flat[i]
    ax.set_yticks([])
    ax.set_xticks([])

plt.tight_layout()
plt.subplots_adjust(wspace=0.2, hspace=0)
'''

# figure2: detection results
plt.rc('font', family='Times New Roman')
plt.rcParams["font.weight"] = 'bold'
plt.rcParams["axes.labelweight"] = 'bold'
plt.rcParams["mathtext.fontset"] = 'stix'

Y, X = img_input.shape

font1 = {'family': 'Times New Roman',
         'weight': 'bold',
         'size': 24,
         }

fig, ax = plt.subplots(1, 1, figsize=(12, 5))
ax.tick_params(length=4, width=1, labelsize=20)
imshow(ax, img_input)
ax.contour(img_border, colors='darkorange')
ax.set_xticks(np.linspace(0, X, 5))
ax.set_xticklabels(['0', '90', '180', '270', '360'])
ax.set_yticks(np.linspace(0, Y, 3))
ax.set_yticklabels(['sin(-90$^\circ$)', '0', 'sin(90$^\circ$)'], fontsize=16)

ax.set_xlabel("Longitude (degree)", fontsize=24)
lat = location[4, :]
lon = location[5, :]
lat1 = (np.sin(lat/180*np.pi) + np.sin(np.pi/3))/(np.sin(np.pi/3)*2)*Y
lon1 = lon*10

for i in range(num):

    ax.text(lon1[i]+80, lat1[i], i+1, fontsize=18)

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)

t_e = time.time()
elapsed_time = t_e - t_s
print(f"runtime: {elapsed_time:.4f} second")
