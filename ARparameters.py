# -*- coding: utf-8 -*-
"""
Getting the parameters of active region

"""
import numpy as np
from skimage import measure
from math import erf
import cv2 as cv
# ---------------------------------------------------------------------------------
rs = 695.5  # 单位：Mm
microH = 2*np.pi*rs**2 / (10**6)  # micro hemisphere
# hign resolution synoptic magnetogram
dx = 360/3600  # degree/pixel,longitude for each pixel
dy = 2/1440  # sin-latitude for each pixel

global pixel_area_Mm, pixel_area_mH
pixel_area_Mm = (dx/180*np.pi)*dy*rs**2  # unit：Mm^2
# unit：micro of the solar hemisphere area
pixel_area_mH = pixel_area_Mm / microH


def GetCentroid(AR, locat):
    """
    Parameters
    ----------
    AR : array
        values of the all pixels of a AR (or its one polarity)
    locat : array
        coordinates of all pixels

    Returns
    -------
    flux weighted center
    lat : latitude
    lon : longitude

    """
    rowN = 1248  # row number of the synoptic map after removing polar field
    colN = 3600  # column number of the map

    # flux weighted centroid
    row = np.sum(AR / np.sum(AR) * locat[0])
    col = np.sum(AR / np.sum(AR) * locat[1])

    # range of map without polar field: -60 60 deg
    sincita = row / rowN * 3**0.5 - 3**0.5/2
    lat = np.arcsin(sincita) * 180 / np.pi
    lon = col / colN * 360

    return lat, lon


def ARLocat(img_label, img_input):
    '''
    get location of all ARs: 
        AR positive polarity location:rows 0-1
        AR negative polarity location:rows 2-3
        AR location:rows 4-5

    '''
    location = np.zeros((6, np.max(img_label)))
    for i in range(1, np.max(img_label)+1):
        locat = np.where(img_label == i)
        AR = img_input[locat]
        pos = AR[np.where(AR > 0)]
        neg = AR[np.where(AR < 0)]

        # of whole AR
        location[4, i-1], location[5, i-1] = GetCentroid(np.abs(AR), locat)

        locat_pos = np.where((img_input > 0) * (img_label == i))
        locat_neg = np.where((img_input < 0) * (img_label == i))

        location[0, i-1], location[1, i-1] = GetCentroid(pos, locat_pos)
        location[2, i-1], location[3, i-1] = GetCentroid(abs(neg), locat_neg)

    return location


def Distance(loc):
    """
    get the distance between two points

    Parameters
    ----------
    loc : [lat0, lon0, lat1, lon1],
    location of two points, unit：degree（°）

    Returns
    -------
    delta : distance between two points, unit：degree（°）

    """
    # change the unit from degree to radian
    lat0, lon0, lat1, lon1 = loc / 180 * np.pi
    #  formular: D=R*arccos(cos(lat0)*cos(lat1)*cos(lon1-lon0) + sin(lat1)*sin(lat0))
    delta = np.arccos(np.cos(lat0) * np.cos(lat1) * np.cos(lon1 - lon0)
                      + np.sin(lat1)*np.sin(lat0))
    # change the unit from radian to degree
    delta = delta / np.pi * 180

    return delta


def Tilt(loc):
    """
    calculate the tilt angle

    Parameters
    ----------
    loc : [lat0, lon0, lat1, lon1] ,(4,)
    location of two points, unit：degree（°）

    Returns
    -------
    tilt : float, unit：deg（°）

    """
    # np.cos,sin 使用弧度制，输入输出的角度使用°，要变换单位

    # 太阳半径
    # rs = 695.5  #单位：Mm
    lat0, lon0, lat1, lon1, latw, lonw = loc / 180 * np.pi
    # 球面两点距离公式 D=R*arccos(cos(lat0)*cos(lat1)*cos(lon1-lon0)
    #                    +sin(lat1)*sin(lat0))
    tilt = np.arctan(-(lat0 - lat1)/((lon0-lon1)*np.cos(latw)))
    tilt = tilt / np.pi * 180

    return tilt


def ARArea(img_label, img_input):
    '''
    Parameters
    ----------
    img_label : array
        img labeled with measure.label
    input_file : opened fits file
         magnetogram file

    Returns
    -------
    area : array
        two plarity area of each AR in the magnetogram; pos,neg
    '''
    # pixel_area_Mm = 1.172570975148625 Mm^2
    # pixel_area_maH = 0.38580246913580246

    area = np.zeros((2, np.max(img_label)))

    for i in range(np.max(img_label)):
        ar = img_input*(img_label == (i+1))
        arn = ar[ar < 0]
        arp = ar[ar > 0]
        # area of postive polarity
        area[0, i] = len(arp)*pixel_area_mH
        # area of negative polarity
        area[1, i] = len(arn)*pixel_area_mH

    return area


def ARFlux(img_label, img_input):
    '''
    get the flux of two polarities of AR;   unit：Mx

    '''
    flux = np.zeros((2, np.max(img_label)))

    for i in range(1, np.max(img_label)+1):
        AR = img_input[np.where(img_label == i)]
        pos = AR[np.where(AR > 0)]
        neg = AR[np.where(AR < 0)]
        # 1 Gauss*Mm^2 = 10^16 Mx
        # positive polarity
        flux[0, i-1] = np.sum(pos) * pixel_area_Mm * 10**16
        # negative polarity
        flux[1, i-1] = np.sum(neg) * pixel_area_Mm * 10**16

    return flux


def ARBmax(img_label, img_input):
    '''
    get the maximum magnetic field of AR;   unit：G
    '''
    Bmax = np.zeros(np.max(img_label))

    for i in range(1, np.max(img_label)+1):
        AR = img_input[np.where(img_label == i)]
        Bmax[i-1] = np.max(abs(AR))

    return Bmax


def FDF(img_label, img_input, lamr):
    """
    calculate the AR final dipole fields

    Parameters
    ----------
    img_label : array
        label of detected ARs
    img_input : array
        original image 

    Returns
    -------
    final dipole field of all detected ARs

    """

    # AR number
    N = np.max(img_label)

    M1, N1 = np.shape(img_label)

    # parameters
    dx = (360/3600)/180*np.pi
    dy = 2/1440
    A = 0.21
    lamr = lamr/180*np.pi  # change unit from degree to rad

    # calculate the FDF based on Wang 2021
    fdf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]
        ar2 = ar != 0

        # balance the flux
        # Wang 2021, improve the weak polarity
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        if p > n:
            index1 = (img_label == i+1)*(img_input < 0)
            ar[index1] = ar[index1]*p/n
        else:
            index1 = (img_label == i+1)*(img_input > 0)
            ar[index1] = ar[index1]*n/p
        '''
        # yeates 2020
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        a = (p+n)/(2*p)
        b = (p+n)/(2*n)
        ar[ar > 0] = ar[ar > 0]*a
        ar[ar < 0] = ar[ar < 0]*b
        '''
        fdf1 = 0
        for j in range(M1):
            for k in range(N1):
                if ar2[j, k]:
                    lat = np.arcsin((j+0.5)/M1*3**0.5 - 0.5*3**0.5)
                    fdf1 += ar[j, k] * erf(lat/(lamr*2**0.5))
        fdf1 = A*fdf1*dx*dy
        fdf.append(fdf1)

    return np.array(fdf)


def IDF(img_label, img_input):
    """
    calculate the AR initial dipole field

    Parameters
    ----------
    img_label : array
        label of detection results
    img_input : array
        original image 

    Returns
    -------
    final dipole field of all detected ARs

    """

    # AR number
    N = np.max(img_label)

    M1, N1 = np.shape(img_label)

    # parameters
    dy = 2/1440
    dx = (360/3600)/180*np.pi

    idf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]
        ar2 = ar != 0

        # balance the flux
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        if p > n:
            index1 = (img_label == i+1)*(img_input < 0)
            ar[index1] = ar[index1]*p/n
        else:
            index1 = (img_label == i+1)*(img_input > 0)
            ar[index1] = ar[index1]*n/p

        idf1 = 0
        for j in range(M1):
            for k in range(N1):
                # speed up the program
                if ar2[j, k]:
                    sinlat = (j+0.5)/M1*3**0.5 - 0.5*3**0.5
                    idf1 += ar[j, k]*sinlat

        idf1 = 3/(4*np.pi)*idf1*dy*dx
        idf.append(idf1)

    return np.array(idf)


def FDFBMR(flux, loc, lamdr):
    """
    calculate the IDF and FDF with BMR approxiamtion

    Parameters
    ----------
    flux : array of AR flux
    loc : aray of AR location
    lamdr : transport parameter

    Returns
    -------
    idf_bmr,   fdf_bmr

    """
    lat = loc[4, :]
    lat = lat*np.pi/180
    num = flux.shape[1]
    flux2 = np.zeros(num)
    for i in range(num):
        flux2[i] = np.max(abs(flux[:, i]))

    '''
    # polarity separation
    d = Distance(loc[0:4, :])
    d = d*np.pi/180
    # tilt angle
    tilt = abs(Tilt(loc))
    tilt = tilt*np.pi/180
    dlat = d*np.sin(tilt)
    sign = np.sign(loc[0, :]-loc[2, :])
    '''
    dlat = (loc[0, :]-loc[2, :])*np.pi/180

    # radius of the sun
    rs = 695.5
    idf_bmr = 3/(4*np.pi*rs**2)*flux2*dlat*np.cos(lat)*1e-16

    d = Distance(loc[0:4, :])

    lamdr = (lamdr**2+(d/2)**2)**0.5
    lamdr = lamdr/180*np.pi
    a = (2/np.pi)**0.5*(8/9)
    finf = a/lamdr*np.exp(-lat**2/(2*lamdr**2))
    fdf_bmr = finf*idf_bmr

    return idf_bmr, fdf_bmr
