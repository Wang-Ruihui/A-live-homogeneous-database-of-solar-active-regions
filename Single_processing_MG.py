# -*- coding: utf-8 -*-
"""
process the single magnetogram
"""
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
import scipy
from ARdetection import ARdetection, imshow, Get_ARi
from ARparameters import ARArea, ARFlux, ARLocat, Distance, FDF, IDF, FDFrs, FDFBMR
# ------------------------------------------------------------------------------
cr = 2148
# pathMDI1 = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.'
# #pathMDI1 = 'D:/python program/活动区识别/data/SynopticMap/Repeated Data/MDI/mdi.synoptic_Mr_96m.'
# pathMDI2 = '.data.fits'
# File = pathMDI1 + str(cr) + pathMDI2
# File_last = pathMDI1 + str(cr-1) + pathMDI2
# File_next = pathMDI1 + str(cr+1) + pathMDI2

pathHMI1 = 'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.'
pathHMI2 = '.synopMr.fits'
File = pathHMI1 + str(cr) + pathHMI2
File_last = pathHMI1 + str(cr-1) + pathHMI2
File_next = pathHMI1 + str(cr+1) + pathHMI2

# use the length of file name to judge its source: HMI or MDI
Filename = File.split('/')[-1]
if len(Filename) > 35:
    filetype = 'HMI'
else:
    filetype = 'MDI'

rt1 = 0.85
rt2 = 1


# read fits file and get the synoptic magnetogram data
img_input, module1, module2, module3, module4, module42, module5,\
    result, img_label, img_border = \
    Get_ARi(File, filetype, File_last, File_next, rt1, rt2)

# results:img_input, module1, module2, module3, module4,module5,  result, img_label, img_border

area = ARArea(img_label, img_input)
flux = ARFlux(img_label, img_input)
number = np.max(img_label)

locat = ARLocat(img_label, img_input)

distance = np.zeros(number)
for i in range(number):
    loc = locat[0:4, i]
    distance[i] = Distance(loc)

lamdr = 7.17
fdf = FDF(img_label, img_input, lamdr)
idf = IDF(img_label, img_input)
fdfrs = FDFrs(flux, locat[4, :], 1, lamdr)
idfbmr, fdfbmr = FDFBMR(flux, locat, lamdr)
# ------------------------------------------------------------------------------

plt.rcParams.update({
    'font.family': 'Arial',
    "font.weight": 'normal',
    "axes.labelweight": 'bold',
    'mathtext.fontset': 'stix'
})

Y, X = img_input.shape
# 识别后图像
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
ax.set_yticklabels(['sin(-60$^\circ$)', '0', 'sin(60$^\circ$)'], fontsize=16)

ax.set_xlabel("Longitude (degree)", fontsize=24)
lat = locat[4, :]
lon = locat[5, :]
lat1 = (np.sin(lat/180*np.pi) + np.sin(np.pi/3))/(np.sin(np.pi/3)*2)*Y
lon1 = lon*10

for i in range(number):
    ax.text(lon1[i]+80, lat1[i], i+1, fontsize=18)

# 图像紧凑分布，去除周围空白和图像间隔
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
