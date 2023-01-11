# -*- coding: utf-8 -*-
"""
Getting the parameters of active regions

"""
import numpy as np
# ---------------------------------------------------------------------------------


def GetCentroid(AR, locat):
    """
    get the flux weighted center of any AR or any polarity of an AR
    
    Parameters
    ----------
    AR : all pixels of each AR or each polarity
    locat : the location of each pixels inputted

    Returns
    -------
    flux weighted center
    latitude and longitude

    """
    rowN = 1248  # the number of total rows of synoptic magnetograms
    colN = 3600  # the number of total columns of synoptic magnetograms
    
    # location of flux weighted centroid
    row = np.sum(AR / np.sum(AR) * locat[0])
    col = np.sum(AR / np.sum(AR) * locat[1])

    # change to lat and lon
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

        # location of entire AR
        location[4, i-1], location[5, i-1] = GetCentroid(np.abs(AR), locat)

        # location of positive and negative polarities of AR
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


def ARArea(img_label, img_input):
    '''
    get the area of two polarities of AR
    '''
    # pixel_area_Mm = 1.172570975148625 Mm^2
    # pixel_area_maH = 0.38580246913580246
    
    # solar radius
    rs = 695.5  # unit：Mm
    microH = 2*np.pi*rs**2 / (10**6)  # micro hemisphere of solar surface
    # hign resolution synoptic magnetogram
    dx = 360/3600  # longitude of each pixel
    dy = 2/1440  # sin latitude of each pixel

    global pixel_area_Mm
    pixel_area_Mm = (dx/180*np.pi)*dy*rs**2  # unit：Mm^2
    pixel_area_mH = pixel_area_Mm / microH  # unit：micro hemisphere(mHem)

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
        # flux of positive polarity
        flux[0, i-1] = np.sum(pos) * pixel_area_Mm * 10**16 
        # flux of negative polarity
        flux[1, i-1] = np.sum(neg) * pixel_area_Mm * 10**16

    return flux
