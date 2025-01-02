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
microH = 2*np.pi*rs**2 / (10**6)  # 太阳半球面积的百万分之一,micro hemisphere
# hign resolution synoptic magnetogram
dx = 360/3600  # degree/pixel,综合磁图上每个像素x方向对映的经度
dy = 2/1440  # 综合磁图上每个像素y方向对映的正弦纬度

global pixel_area_Mm, pixel_area_mH
pixel_area_Mm = (dx/180*np.pi)*dy*rs**2  # 单位：Mm^2
pixel_area_mH = pixel_area_Mm / microH  # 单位：太阳半球面积的百万分之一


def GetCentroid(AR, locat):
    """
    Parameters
    ----------
    AR : array
        每个活动区所有像素磁场值，所有正磁场值或所有负磁场值
    locat : array
        每个活动区对应输入的AR类型的像素的坐标

    Returns
    -------
    flux weighted center
    lat : 纬度
    lon : 经度

    """
    rowN = 1248  # 修改后的综合磁图的总行数
    colN = 3600  # 综合磁图总列数

    '''
    #几何中心
    row = np.average(locat[0])
    col = np.average(locat[1])

    '''
    # 重心的矩阵下标,以磁通量为权重
    row = np.sum(AR*locat[0])/np.sum(AR)
    col = np.sum(AR*locat[1])/np.sum(AR)

    # 重心的经纬度
    # 纬度对应的正弦值，没有高纬地区，磁图对应纬度范围+-60
    sincita = (row+0.5) / rowN * 3**0.5 - 3**0.5/2
    lat = np.arcsin(sincita) * 180 / np.pi
    lon = (col+0.5) / colN * 360

    return lat, lon


def ARLocat(img_label, img_input):
    '''
    img_input : array
        输入的原图像
    img_label : array
        活动区标记后的图像

    Returns
    -------
    location : array
        所有活动区的经纬度位置
        每一列表示一个活动区数据；
        6行
        0，1行表示活动区正极区经纬度位置(0纬度，1经度)，
        2和3表示活动区负极区经纬度位置，
        4和5表示活动区整体经纬度位置；
        每一部分都是第一行表示活动区纬度，第二行表示活动区经度

    '''
    location = np.zeros((6, np.max(img_label)))
    for i in range(1, np.max(img_label)+1):
        locat = np.where(img_label == i)
        AR = img_input[locat]
        pos = AR[np.where(AR > 0)]
        neg = AR[np.where(AR < 0)]

        # 活动区整体位置
        location[4, i-1], location[5, i-1] = GetCentroid(np.abs(AR), locat)

        # 获取活动区正负极像素的下标
        '''
        AR2 = np.zeros(img_label.shape)
        AR2[np.where(img_label == i)] = 1
        AR2 = AR2 * img_input
        locat_pos = np.where(AR2 > 0)
        locat_neg = np.where(AR2 < 0)
        '''
        locat_pos = np.where((img_input > 0) * (img_label == i))
        locat_neg = np.where((img_input < 0) * (img_label == i))

        location[0, i-1], location[1, i-1] = GetCentroid(pos, locat_pos)
        location[2, i-1], location[3, i-1] = GetCentroid(abs(neg), locat_neg)

    return location


def Distance(loc):
    """
    求解球面两点的球心角

    Parameters
    ----------
    loc : [lat0, lon0, lat1, lon1] ,(4,)
        两点的纬度、经度坐标 ,单位：度（°）

    Returns
    -------
    delta : float
        球面两点的球心角, 单位：度（°）

    """
    # np.cos,sin 使用弧度制，输入输出的角度使用°，要变换单位

    # 太阳半径
    # rs = 695.5  #单位：Mm
    lat0, lon0, lat1, lon1 = loc / 180 * np.pi
    # 球面两点距离公式 D=R*arccos(cos(lat0)*cos(lat1)*cos(lon1-lon0)
    #                    +sin(lat1)*sin(lat0))
    delta = np.arccos(np.cos(lat0) * np.cos(lat1) * np.cos(lon1 - lon0)
                      + np.sin(lat1)*np.sin(lat0))
    delta = delta / np.pi * 180

    return delta


def Tilt(loc):
    """
    求解活动区两级倾角

    Parameters
    ----------
    loc : [lat0, lon0, lat1, lon1] ,(4,)
        两极的纬度、经度坐标 ,单位：度（°）

    Returns
    -------
    tilt : float
        球面两点的倾角, 单位：度（°）

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
    # 面积计算，算算
    """
    properties = measure.regionprops(img_label)
    for i in range(len(properties)):
        prop = properties[i]
        area[i] = prop.area * pixel_area_mH
    """

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
    Parameters
    ----------
    img_label : array
        活动区标记后的图像矩阵
    img_input : array
        处理后的图像输入矩阵，1248*3600

    Returns
    -------
    flux : array
        磁通量矩阵，第一行为正磁通量，第二行为负磁通量；每一列为一个活动区

    '''
    # 单位：Mx 麦克斯韦
    flux = np.zeros((2, np.max(img_label)))

    for i in range(1, np.max(img_label)+1):
        AR = img_input[np.where(img_label == i)]
        pos = AR[np.where(AR > 0)]
        neg = AR[np.where(AR < 0)]
        # 1 Gauss*Mm^2 = 10^16 Mx
        flux[0, i-1] = np.sum(pos) * pixel_area_Mm * 10**16  # 正磁通量
        flux[1, i-1] = np.sum(neg) * pixel_area_Mm * 10**16  # 负磁通量

    return flux


def ARBmax(img_label, img_input):
    '''
    get the flux of two polarities of AR;   unit：Mx
    '''
    Bmax = np.zeros(np.max(img_label))

    for i in range(1, np.max(img_label)+1):
        AR = img_input[np.where(img_label == i)]
        Bmax[i-1] = np.max(abs(AR))

    return Bmax


def FDF(img_label, img_input, lamr):
    """
    calculate the final dipole magnetic field of ARs

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
    dx = (360/3600)/180*np.pi
    dy = 2/1440
    A = 0.21
    lamr = lamr/180*np.pi  # change unit from degree to rad

    # calculate the FDF based on Wang, 2021
    fdf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]

        # balance the flux
        # Wang 2021
        # pos = ar[ar > 0]
        # neg = ar[ar < 0]
        # p = np.sum(pos)
        # n = np.sum(abs(neg))
        # if p > n:
        #     index1 = (img_label == i+1)*(img_input < 0)
        #     ar[index1] = ar[index1]*p/n
        # else:
        #     index1 = (img_label == i+1)*(img_input > 0)
        #     ar[index1] = ar[index1]*n/p

        # yeates 2020
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        a = (p+n)/(2*p)
        b = (p+n)/(2*n)
        ar[ar > 0] = ar[ar > 0]*a
        ar[ar < 0] = ar[ar < 0]*b

        fdf1 = 0
        for j in range(M1):
            lat = np.arcsin((j+0.5)/M1*3**0.5 - 0.5*3**0.5)
            fdf1 += np.sum(ar[j, :]) * erf(lat/(lamr*2**0.5))
        fdf1 = A*fdf1*dx*dy
        fdf.append(fdf1)

    return np.array(fdf)


def IDF(img_label, img_input):
    """
    calculate the final dipole magnetic field of ARs

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

    # calculate the FDF based on Wang, 2021
    idf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]

        # balance the flux
        # # wang 2020
        # pos = ar[ar > 0]
        # neg = ar[ar < 0]
        # p = np.sum(pos)
        # n = np.sum(abs(neg))
        # if p > n:
        #     index1 = (img_label == i+1)*(img_input < 0)
        #     ar[index1] = ar[index1]*p/n
        # else:
        #     index1 = (img_label == i+1)*(img_input > 0)
        #     ar[index1] = ar[index1]*n/p

        # yeates 2020
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        a = (p+n)/(2*p)
        b = (p+n)/(2*n)
        ar[ar > 0] = ar[ar > 0]*a
        ar[ar < 0] = ar[ar < 0]*b

        idf1 = 0
        for j in range(M1):
            sinlat = (j+0.5)/M1*3**0.5 - 0.5*3**0.5
            idf1 += np.sum(ar[j, :])*sinlat

        idf1 = 3/(4*np.pi)*idf1*dy*dx
        idf.append(idf1)

    return np.array(idf)


def FDFsm(img_label, img_input, lamr):
    """
    calculate the final dipole magnetic field of ARs

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

    # size of small picture
    M2 = 156
    N2 = 360
    # parameters
    dx = (360/360)/180*np.pi
    dy = 2/180
    A = 0.21
    lamr = lamr/180*np.pi  # change unit from degree to rad

    # calculate the FDF based on Wang, 2021
    fdf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]

        ar_sml = cv.resize(ar, (N2, M2), interpolation=cv.INTER_LINEAR)
        # # balance the flux
        # # wang 2021
        # pos = ar_sml[ar_sml > 0]
        # neg = ar_sml[ar_sml < 0]
        # p = np.sum(pos)
        # n = np.sum(abs(neg))
        # if p > n:
        #     index1 = ar_sml < 0
        #     ar_sml[index1] = ar_sml[index1]*p/n
        # else:
        #     index1 = ar_sml > 0
        #     ar_sml[index1] = ar_sml[index1]*n/p

        # yeates 2020
        indexp = ar_sml > 0
        indexn = ar_sml < 0
        pos = ar_sml[indexp]
        neg = ar_sml[indexn]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        a = (p+n)/(2*p)
        b = (p+n)/(2*n)
        ar_sml[indexp] = ar_sml[indexp]*a
        ar_sml[indexn] = ar_sml[indexn]*b

        fdf1 = 0
        for j in range(M2):
            lat = np.arcsin((j+0.5)/M2*3**0.5 - 0.5*3**0.5)
            fdf1 += np.sum(ar_sml[j, :]) * erf(lat/(lamr*2**0.5))
        fdf1 = A*fdf1*dx*dy
        fdf.append(fdf1)

    return np.array(fdf)


def IDFsm(img_label, img_input):
    """
    calculate the final dipole magnetic field of ARs

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

    # size of small picture
    M2 = 156
    N2 = 360
    # parameters
    dx = (360/360)/180*np.pi
    dy = 2/180

    # calculate the FDF based on Wang, 2021
    idf = []
    for i in range(N):
        ar = np.zeros((M1, N1))
        index = (img_label == i+1)
        ar[index] = img_input[index]

        # reduce the size
        ar_sml = cv.resize(ar, (N2, M2), interpolation=cv.INTER_LINEAR)

        # # balance the flux
        # wang 2020
        # pos = ar_sml[ar_sml > 0]
        # neg = ar_sml[ar_sml < 0]
        # p = np.sum(pos)
        # n = np.sum(abs(neg))
        # if p > n:
        #     index1 = ar_sml < 0
        #     ar_sml[index1] = ar_sml[index1]*p/n
        # else:
        #     index1 = ar_sml > 0
        #     ar_sml[index1] = ar_sml[index1]*n/p

        # yeates 2020
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p = np.sum(pos)
        n = np.sum(abs(neg))
        a = (p+n)/(2*p)
        b = (p+n)/(2*n)
        ar[ar > 0] = ar[ar > 0]*a
        ar[ar < 0] = ar[ar < 0]*b

        idf1 = 0
        for j in range(M2):
            sinlat = (j+0.5)/M2*3**0.5 - 0.5*3**0.5
            idf1 += np.sum(ar_sml[j, :])*sinlat

        idf1 = 3/(4*np.pi)*idf1*dy*dx
        idf.append(idf1)

    return np.array(idf)


def FDFrs(flux, lat, sign, lamdr):
    # sign is the sign of flux in the northern polarity
    lat = lat*np.pi/180
    num = flux.shape[1]
    flux2 = np.zeros(num)
    for i in range(num):
        flux2[i] = np.max(abs(flux[:, i]))

    # 极间距
    # Lemerle 2015
    #d = np.exp((0.46 + 0.42*(np.log(flux2)/np.log(10)-21))*np.log(10))
    # ours
    d = (flux2/(5.7e20))**(1/1.73)  # flux(d)
    # d = (flux2/(np.exp(47.77)))**(1/1.8)  # flux(d) Cy 23
    # d = (flux2/(np.exp(47.52)))**(1/1.83)  # flux(d) smooth map
    # d = np.exp(0.305*np.log(flux2)-13.83)  # d(flux)
    # yeates2020
    #d = (flux2/(3.2e20))**(1/1.8)
    d = d*np.pi/180
    # 倾角
    #tilt = 0.5*abs(lat)
    #tilt = abs(0.42*lat+0.37)
    tilt = abs(0.42*lat)  # cycle23, all AR,no b
    # tilt = abs(0.41*lat)  # smooth map
    dlat = d*np.sin(tilt)
    # radius of the sun
    rs = 695.5
    idf_rs = sign*3/(4*np.pi*rs**2)*flux2*dlat*np.cos(lat)*1e-16

    lamdr = lamdr/180*np.pi
    a = (2/np.pi)**0.5*(8/9)
    finf = a/lamdr*np.exp(-lat**2/(2*lamdr**2))
    fdf_rs = finf*idf_rs

    return fdf_rs


def FDFBMR(flux, loc, lamdr):
    # sign is the sign of flux in the northern polarity
    lat = loc[4, :]
    lat = lat*np.pi/180
    num = flux.shape[1]
    flux2 = np.zeros(num)
    # # flux balance wang 2020 (increase the weak polarity)
    # for i in range(num):
    #     flux2[i] = np.max(abs(flux[:, i]))

    # flux balance Yeates 2020 (keep the total unsigned flux)
    for i in range(num):
        flux2[i] = np.mean(abs(flux[:, i]))
    '''
    # 极间距
    d = Distance(loc[0:4, :])
    d = d*np.pi/180
    # 倾角
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
