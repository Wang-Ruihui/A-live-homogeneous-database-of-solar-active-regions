# -*- coding: utf-8 -*-
"""
batch process the magnetogram

生成活动区参数矩阵，保存到文件File2
"""

import numpy as np
from astropy.io import fits
from skimage import measure
import os
from ARdetection import Get_ARP, Get_ARi

# ------------------------------------------------------------------------------


def data_missing(file):
    """
    Determine whether the file is missing data
    if there is nan around equator of the map, the file is missing data

    Parameters
    ----------
    file : filepath of fits file

    Returns
    -------
    bool
        True for missing，False for no missing
    """
    input_file = fits.open(file)
    input_img = input_file[0].data
    input_file.close()

    num_nan = len(np.where(np.isnan(input_img[720, :]))[0])

    if num_nan >= 1:
        flag = True
    else:
        flag = False

    return flag


# --------------------------------------------------------------------------------
if __name__ == '__main__':

    data_missing_file = []

    # batch process the 23 and 24 solar cycle magneticgram
    Filepath = 'D:\\python program\\活动区识别\\SynopMr'
    lists1 = os.listdir(Filepath)

    AR_parameters = np.zeros((17, 1))
    file_last = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'
    cr = 1909
    rt1 = 0.85
    rt2 = 1
    lamr = 7.17

    for j in range(len(lists1)):
        path = lists1[j]
        filepath_new = Filepath+'\\'+path  # folder of each instrument
        lists2 = os.listdir(filepath_new)

        fileNum = len(lists2)
        print(fileNum)

        for i in range(fileNum):
            filename1 = lists2[i]
            file = os.path.join(filepath_new, filename1)  # file of each map
            if i+1 < fileNum:
                file_next = os.path.join(filepath_new, lists2[i+1])
            elif j == 0:
                file_next = 'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.2097.synopMr.fits'
            else:
                file_next = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'
            # determine whether the file is missing data. if true, label
            if data_missing(file):
                data_missing_file.append(filename1)

            # determine the file type, MDI or HMI
            if len(filename1) > 35:
                filetype = 'HMI'
            else:
                filetype = 'MDI'

            parameters = Get_ARP(file, filetype, file_last,
                                 file_next, rt1, rt2, lamr)
            AR_parameters = np.concatenate((AR_parameters, parameters), axis=1)

            cr += 1
            if cr == 1938:
                cr = 1941

            file_last = file

            print(filename1+' done')
            print('      ')

        print(path.split()[1] + ' magnetogram done')

    # save the data
    AR_parameters = np.delete(AR_parameters, 0, axis=1)  # delete the first col
    File2 = 'D:/python program/活动区识别/result data/Remove12/AR_Parameters_v1(lamr7.17,BL_yeates).npy'
    np.save(File2, AR_parameters)
