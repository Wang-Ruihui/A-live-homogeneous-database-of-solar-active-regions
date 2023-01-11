# -*- coding: utf-8 -*-
"""
Automatic detection of solar active regions from synoptic magnetograms

"""

import numpy as np
import cv2
from skimage import measure, morphology

# ------------------------------------------------------------------------------


def imshow(ax, data):
    """
    show image in 3 cigma way
    """
    mean = np.average(data)
    std = np.std(data)
    Y, X = data.shape
    vmin = mean-3*std
    vmax = mean+3*std
    if vmax > np.max(data):
        vmax = np.max(data)
    if vmin < np.min(data):
        vmin = np.min(data)

    # cmap='gray'
    ax.imshow(data, cmap='gray', vmin=vmin, vmax=vmax,
              origin='lower')  


def adaptiveThreshold(img_input, size, const):
    """
    adaptive threshold segmentation

    """
    img_smt = cv2.GaussianBlur(img_input, (size, size), 0, borderType=2)
    img_threshold = img_smt + const
    img_result = img_input - img_threshold

    img_result[np.where(img_result > 0)] = 1
    img_result[np.where(img_result < 0)] = 0

    return img_result


def threshold(img_input):
    """
    intensity threshold segmentation
    """

    # img_input data type: 'int16'
    img_input1 = np.abs(img_input).astype("int16")

    # adaptive threshold
    img_thresh = adaptiveThreshold(img_input1, 501, 10)

    return img_thresh

# ------------------------------------------------------------------------------


class Point(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def getX(self):
        return self.x

    def getY(self):
        return self.y


def selectConnects(p):
    # connection types,p=0: 4 connection; p=1：8 connection;
    if p != 0:
        connects = [Point(-1, -1), Point(0, -1), Point(1, -1), Point(1, 0),
                    Point(1, 1), Point(0, 1), Point(-1, 1), Point(-1, 0)]
    else:
        connects = [Point(0, -1), Point(1, 0), Point(0, 1), Point(-1, 0)]
    return connects


def findseeds(array, origion_img, thresh):
    # find seeds used in the region growing. Each seed should be stronger than 50G
    seedlist = []
    # the coordinates of remain pixels in array
    index = np.where(array == 1)
    # the values of remain pixels in array
    value = origion_img[array == 1]

    # get the coordinates of pixels stronger than 50G
    seedindex = np.where(abs(value) > 50)[0]
    for i in range(len(seedindex)):
        index1 = seedindex[i]
        point = Point(index[0][index1], index[1][index1])
        seedlist.append(point)

    return seedlist


def regionGrow(img, seeds, thresh, p=1):
    row, col = img.shape
    # image after region growing
    seedMark = np.zeros(img.shape, 'uint8')
    seedList = seeds[:]
    label = 1
    # default p=1
    connects = selectConnects(p)
    while (len(seedList) > 0):
        currentPoint = seedList.pop(0)

        seedMark[currentPoint.x, currentPoint.y] = label
        if p == 1:
            n = 8
        else:
            n = 4

        for i in range(n):
            # the coordinates of neighbour pixels
            tmpX = currentPoint.x + connects[i].x
            tmpY = currentPoint.y + connects[i].y

            # if the neignbour pixel are beyond the image, skip it
            if tmpX < 0 or tmpY < 0 or tmpX >= row or tmpY >= col:
                continue
            # if the neignbour pixel are not labelled, judge its value
            elif seedMark[tmpX, tmpY] == 0:
                Bval = abs(img[tmpX, tmpY])
                if Bval > thresh:
                    seedMark[tmpX, tmpY] = label
                    seedList.append(Point(tmpX, tmpY))
    return seedMark

# ------------------------------------------------------------------------------


def RemoveUnipolar(img_input, img_label):
    # remove unipolar regions
    # relabel the remain ARs
    flag = 1
    for i in range(1, np.max(img_label) + 1):
        AR = img_input[np.where(img_label == i)]
        pos = AR[np.where(AR > 0)]
        neg = AR[np.where(AR < 0)]
        # totoal flux
        FluxSign = np.sum(pos) + np.sum(neg)
        # total ungisn flux
        FluxUnSign = np.sum(pos) + np.abs(np.sum(neg))

        if abs(FluxSign) >= 0.5 * FluxUnSign:
            img_label[np.where(img_label == i)] = 0
        else:
            img_label[np.where(img_label == i)] = flag
            flag = flag + 1

    return img_label

# ------------------------------------------------------------------------------


def ARdetection(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4):
    """
    img_input: input FITS file
    Kernel1  : closing operation kernel in module2
    Kernel2  : opening operation kernel in module2
    Thresh   : region grow threshold
    Kernel3  : closing operation kernel in module4
    Size     : the smallest size of AR; any region smaller than it will be removed.
    Kernel4  : dilating operation kernel in module5

    return:  result after each module, detection result, 
            boundary and labels of all ARs
    """
    # module1: intensity threshold segmentation, adaptive threshold
    module1 = threshold(img_input)

    # module2: morphological closing operation and opening operation to
    #          remove small magnetic segments and get the kernel pixels of ARs
    module2 = cv2.morphologyEx(module1, cv2.MORPH_CLOSE, Kernel1)
    module2 = cv2.morphologyEx(module2, cv2.MORPH_OPEN, Kernel2)

    # module3: region growing
    seeds = findseeds(module2, img_input, Thresh)
    module3 = regionGrow(img_input, seeds, Thresh)

    # module4: closing operation and remove small decayed ARs segments
    module4 = cv2.morphologyEx(module3, cv2.MORPH_CLOSE, Kernel3)
    module4 = morphology.remove_small_objects(module4.astype('bool'), Size)

    # module5: merging neighbor regions and removing unipolar regions
    # merge neighbor regions
    module4 = module4.astype("uint8")
    img_shape = cv2.dilate(module4, Kernel4)
    img_label = measure.label(img_shape)
    img_label = img_label * module4

    # remove unipolar regions
    img_label = RemoveUnipolar(img_input, img_label)
    module5 = np.where(img_label >= 1, 1, 0)

    # get the border of each merged AR
    img_shape1 = np.where(img_label >= 1, 1, 0).astype('uint8')
    img_shape2 = cv2.dilate(img_shape1, Kernel4)
    img_border = img_shape2 - cv2.erode(img_shape2, np.ones((3, 3)))

    # final detection result
    result = module5 * img_input

    return module1, module2, module3, module4, module5, result, \
        img_label, img_border
