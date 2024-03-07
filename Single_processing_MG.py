# -*- coding: utf-8 -*-
"""
process the single magnetogram
"""
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
import scipy
from ARdetection import ARdetection, imshow, Get_ARi
from ARparameters import ARArea, ARFlux, ARLocat, Distance, FDF, IDF, FDFBMR
# ------------------------------------------------------------------------------
# File: image needs processing; File_old:image from the last CR
File = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.1976.data.fits'
File_old = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.1975.data.fits'

# use the length of file name to judge its source: HMI or MDI
Filename = File.split('/')[-1]
if len(Filename) > 35:
    filetype = 'HMI'
else:
    filetype = 'MDI'

rt = 0.7

img_input, module1, module2, module3, module4, module42, module5,  result, \
    img_label, img_border = Get_ARi(File, filetype, File_old, rt)


area = ARArea(img_label, img_input)
flux = ARFlux(img_label, img_input)
number = np.max(img_label)

locat = ARLocat(img_label, img_input)

distance = np.zeros(number)
for i in range(number):
    loc = locat[0:4, i]
    distance[i] = Distance(loc)

lamdr = 5
fdf = FDF(img_label, img_input, lamdr)
idf = IDF(img_label, img_input)
fdfrs = FDFrs(flux, locat[4, :], 1, lamdr)
idfbmr, fdfbmr = FDFBMR(flux, locat, lamdr)
# ------------------------------------------------------------------------------
#
plt.rc('font', family='Times New Roman')
plt.rcParams["font.weight"] = 'bold'
plt.rcParams["axes.labelweight"] = 'bold'
plt.rcParams["mathtext.fontset"] = 'stix'

Y, X = img_input.shape
# 识别后图像
font1 = {'family': 'Times New Roman',
         'weight': 'bold',
         'size': 24,
         }

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
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

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
