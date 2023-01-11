# -*- coding: utf-8 -*-
"""
process the single magnetogram
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from ARdetection import ARdetection, imshow
from ARparameters import ARArea, ARFlux, ARLocat, Distance
# ------------------------------------------------------------------------------
# CR1968, CR2155
File = 'G:/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.1968.data.fits'
#File = 'G:/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.2155.synopMr.fits'

# read fits file and get the synoptic magnetogram data
input_file = fits.open(File)
img_input = input_file[0].data
input_file.close()

# replace the nan in img_input with 0
img_input[np.where(np.isnan(img_input))] = 0

# limit the latitude to (-60,60) or (-pi/3,pi/3) to remove polar magnetic field
rowN = img_input.shape[0]
PR = int((1 - np.sin(np.pi/3))*rowN/2)
img_input = img_input[PR:rowN-PR, :]


# use the length of file name to judge its source: HMI or MDI
Filename = File.split('/')[-1]
if len(Filename) < 38:
    # MDI
    # closing operation kernel in module2
    Kernel1 = np.ones((3, 3))
    # opening operation kernel in module2
    Kernel2 = np.ones((11, 11))
    # region grow threshold
    Thresh = 50
    # closing operation kernel in module4
    Kernel3 = np.ones((5, 5))
    # ARs area threshold in module4
    Size = 351
    # dilating operation kernel in module5
    Kernel4 = np.ones((23, 23))
else:
    # HMI
    Kernel1 = np.ones((3, 3))
    Kernel2 = np.ones((9, 9))
    Thresh = 30
    Kernel3 = np.ones((5, 5))
    Size = 351
    Kernel4 = np.ones((23, 23))

# detect ARs 
module1, module2, module3, module4, module5,  result, img_label, img_border,\
    = ARdetection(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4)


area = ARArea(img_label, img_input)
flux = ARFlux(img_label, img_input)
number = np.max(img_label)

R = []
R.append(number)
R.append(np.sum(area))
R.append(np.sum(abs(flux)))

locat = ARLocat(img_label, img_input)

distance = np.zeros(number)
for i in range(number):
    loc = locat[0:4, i]
    distance[i] = Distance(loc)
# ------------------------------------------------------------------------------
plt.rcParams["axes.labelweight"] = 'bold'
plt.rcParams["mathtext.fontset"] = 'stix'

Y, X = img_input.shape

fig, axes = plt.subplots(7, 1,  figsize=(16, 8))
img = [img_input, module1, module2,
       module3, module4, module5, result]

for i in range(len(axes.flat)):
    ax = axes.flat[i]
    imshow(ax, img[i])
    if i == 0:
        ax.contour(img_border, colors='darkorange')

axes[0].tick_params(length=4, width=1, labelsize=18)
axes[0].set_xticks([])
axes[0].set_yticks(np.linspace(0, Y, 3))
axes[0].set_yticklabels(['sin(-60)', '0', 'sin(60)'])
axes[6].tick_params(length=4, width=1, labelsize=18)
axes[6].set_xticks(np.linspace(0, X, 5))
axes[6].set_xticklabels(['0', '90', '180', '270', '360'])
axes[6].set_yticks(np.linspace(0, Y, 3))
axes[6].set_yticklabels(['sin(-60)', '0', 'sin(60)'])

# remove the ticks of middle panel
for i in range(1, len(axes.flat)-1):
    ax = axes.flat[i]
    ax.set_yticks([])
    ax.set_xticks([])

plt.tight_layout()
plt.subplots_adjust(wspace=0.2, hspace=0)


