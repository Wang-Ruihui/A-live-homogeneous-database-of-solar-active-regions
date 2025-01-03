# -*- coding: utf-8 -*-
"""
Created on Fri May 24 20:57:26 2024

output the detected ARs in the given size

@author: Ruihui Wang
"""
import numpy as np
import cv2 as cv
import os
from sunpy.coordinates.sun import carrington_rotation_time
from ARdetection import Get_ARi
from ARparameters import ARLocat


def FluxBalanceL(ar):
    # balance the flux of ar in lat and lon grid; can not sum all values directly

    # yeates 2020 method
    m, n = np.shape(ar)
    pos = np.zeros((m, n))
    neg = np.zeros((m, n))
    pos[ar > 0] = ar[ar > 0]
    neg[ar < 0] = abs(ar[ar < 0])

    p = 0
    n = 0
    for i in range(m):
        p += np.sum(pos[i, :]*np.sin((i+0.5)*np.pi/180))
        n += np.sum(neg[i, :]*np.sin((i+0.5)*np.pi/180))

    a = (p+n)/(2*p)
    b = (p+n)/(2*n)
    ar[ar > 0] = ar[ar > 0]*a
    ar[ar < 0] = ar[ar < 0]*b

    return ar


def SinlatToLat(imgSinlat):
    # Coordinate transformation: from sinlat to lat

    M, N = imgSinlat.shape
    imgLat = np.zeros((M, N))
    # iteration of equal latitude array
    for i in range(M):
        lat = ((i+0.5)/(M)*180-90)/180*np.pi
        indexlat = (np.sin(lat)+1)/2*(M)

        # linear interpolation
        b = int(indexlat)
        c = indexlat-b
        if b < M-1:
            imgLat[i, :] = (1-c)*imgSinlat[b, :] + c*imgSinlat[b+1, :]
        else:
            imgLat[i, :] = imgSinlat[b, :]

    return imgLat


def OutputARs(directory_path, img_org, img_label, size, CR, CR0, arlat=None, latrange=[-90, 90]):
    """
    output the magnetograms of detected ARs in the given size

    Parameters
    ----------
    img_org : array
        original magnetograms
    img_label : array
        detected AR label of img
    sizes:array,eg. (360,180)
        size of output image
    CR: the CR of map, used to give the time of AR
    CR0: the initial time to give the time of AR
    arlat: The lat list of ARs
    latrange: the lat range of AR that want to output

    Returns
    -------
    None.

    """
    # AR number
    nar = np.max(img_label)
    t0 = carrington_rotation_time(CR0)
    t1 = carrington_rotation_time(CR)
    t2 = carrington_rotation_time(CR+1)

    row = 1440
    col = 3600
    index_lat60 = int((1 - np.sin(np.pi/3))*row/2)

    # list of ARs appearing in the same day
    AR_sameday = []
    # list of ARs appearing on the border of the map
    AR_border = []
    for i in range(nar):
        ar = np.zeros(np.shape(img_org))
        index = (img_label == i+1)
        ar[index] = img_org[index]

        # calculate the date of the AR
        # first calculate the lon
        B1 = abs(ar[index])
        loc = np.where(index)
        col1 = np.sum(B1 * loc[1]) / np.sum(B1)
        lon = (col1+0.5)*0.1

        if arlat is not None:
            if arlat[i] < latrange[0] or arlat[i] > latrange[1]:
                continue

        # set a stricter flux balance limitation for ARs on the border
        if lon > 350 or lon < 10:
            pos = ar[ar > 0]
            neg = ar[ar < 0]
            r = np.sum(pos)/np.sum(abs(neg))
            if r > 2 or r < 0.5:
                AR_border.append(CR)
                AR_border.append(i+1)
                AR_border.append(r)
                continue

        t = t2+lon/360*(t1-t2)
        DeltaDay = (t-t0).value
        DeltaDay = round(DeltaDay)
        # change DeltaDay to str with length 5
        DeltaDay = f"{DeltaDay:05d}"

        # smooth the image before resize to a smaller size
        arSmt = cv.GaussianBlur(ar, (11, 11), 4)

        # restore the latitude range to +-90
        arSmt2 = np.zeros((row, col))
        arSmt2[index_lat60:(row-index_lat60)] = arSmt
        # Coordinate transformation: from sinlat to lat
        arSmt_lat = SinlatToLat(arSmt2)

        arResized = cv.resize(arSmt_lat, size, interpolation=cv.INTER_LINEAR)

        # balance the flux
        arResized = FluxBalanceL(arResized)

        fname = 'day_' + DeltaDay + '.txt'
        file = directory_path + fname
        if os.path.exists(file):
            arResized = arResized + np.loadtxt(file)
            AR_sameday.append(CR)
            AR_sameday.append(i+1)
            AR_sameday.append(DeltaDay)

        np.savetxt(file, arResized)

    return AR_sameday, AR_border


if __name__ == '__main__':
    pathMDI1 = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.'
    pathMDI2 = '.data.fits'
    pathHMI1 = 'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.'
    pathHMI2 = '.synopMr.fits'

    rt1 = 0.85
    rt2 = 1

    # size of output AR maps
    size = (360, 180)

    # set the range of lat, limit to the southern hemisphere
    latrange = [-90, 0]

    ar_sameday = []
    ar_border = []
    directory_path = 'D:\\test\\'

    #File_next = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'
    for i in range(2141, 2160+1):
        cr = i
        # File = pathMDI1 + str(cr) + pathMDI2
        # File_last = pathMDI1 + str(cr-1) + pathMDI2
        # File_next = pathMDI1 + str(cr+1) + pathMDI2

        File = pathHMI1 + str(cr) + pathHMI2
        File_last = pathHMI1 + str(cr-1) + pathHMI2
        File_next = pathHMI1 + str(cr+1) + pathHMI2

        results = Get_ARi(File, 'HMI', File_last, File_next, rt1, rt2)

        img_org = results[0]
        img_label = results[8]
        arlat = ARLocat(img_label, img_org)[4, :]

        # save AR maps
        # take the CR1911 as initial contion, then C1912 is the initial time
        ar_sameday2, ar_border2 = OutputARs(
            directory_path, img_org, img_label, size, cr, 1912, arlat, latrange)
        if len(ar_sameday2) > 0:
            ar_sameday.append(ar_sameday2)
        if len(ar_border2) > 0:
            ar_border.append(ar_border2)
        print(cr)
