# -*- coding: utf-8 -*-
"""
batch process the magnetogram

"""

import numpy as np
from astropy.io import fits
from skimage import measure
import os
from ARdetection import Get_ARP, Get_ARi

# ------------------------------------------------------------------------------


def data_missing(file):
    """
    determine thether there is missing value in the file.
    if there is nan around the equator in the file, then we think is data missing.

    Parameters
    ----------
    file : file path

    Returns
    -------
    bool
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
    # detection of all MDI and HMI maps
    data_missing_file = []

    # batch process the 23 and 24 solar cycle magneticgram
    Filepath = 'D:\\python program\\活动区识别\\SynopMr'
    lists1 = os.listdir(Filepath)

    AR_parameters = np.zeros((17, 1))
    # CR 2077 map used as a map of no ARs
    file_old = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'
    cr = 1909
    rt = 0.7
    lamr = 7

    for j in range(len(lists1)):
        path = lists1[j]
        filepath_new = Filepath+'\\'+path  # path of the fold of each instrument map
        lists2 = os.listdir(filepath_new)

        fileNum = len(lists2)
        print(fileNum)

        for i in range(fileNum):
            filename1 = lists2[i]
            file = os.path.join(filepath_new, filename1)  # path of each map

            # determine whether the file miss data. if yes, then label it
            if data_missing(file):
                data_missing_file.append(filename1)

            # determine the file type，MDI or HMI
            if len(filename1) > 35:
                filetype = 'HMI'
            else:
                filetype = 'MDI'

            parameters = Get_ARP(file, filetype, file_old, lamr, rt)
            AR_parameters = np.concatenate((AR_parameters, parameters), axis=1)

            cr += 1
            if cr == 1938:
                cr = 1941

            file_old = file

            print(filename1+' done')
            print('      ')

        print(path.split()[1] + ' magnetogram done')

    AR_parameters = np.delete(AR_parameters, 0, axis=1)
    File2 = 'D:/python program/活动区识别/result data/onlyclose,rt0.7/Di,Dfbmr/AR Parameters(lamr7.0).npy'
    np.save(File2, AR_parameters)
