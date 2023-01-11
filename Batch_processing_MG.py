# -*- coding: utf-8 -*-
"""
batch process the magnetograms

"""

import numpy as np
from astropy.io import fits
import os
from ARdetection import ARdetection
from ARparameters import ARArea, ARFlux, ARLocat

# ------------------------------------------------------------------------------


def data_missing(file):
    """
    determine whether the data of file is missed
    (if there is any nan value in latitude then the data is missed)

    Parameters
    ----------
    file : fits file

    Returns
    --------
    True for missing; False for complete
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


def GetARP(file, filetype):
    """
    get AR parameters
    file : fits file for magnetogram
    filetype: MDI or HMI

    Returns：param (11*num)
           AR parameters:
           param[0,:]: CR
           param[1,:]: label of each AR
           param[2:4,:]: AR positive polarity location
           param[4:6,:]: AR negative polarity location
           param[6:8,:]: AR location
           param[8:10,:]: AR area
           param[10,:]: AR positive flux
           param[11,:]: AR negative flux
           In AR location, the first row is latitude and the second is longitude.
           Every column is a detected AR.
    """
    input_file = fits.open(file)
    img_input = input_file[0].data
    CR = input_file[0].header['CAR_ROT']
    input_file.close()

    # replace the nan with 0
    img_input[np.where(np.isnan(img_input) == True)] = 0

    # remove polar region and limit the latitude in (-60,60),pi/3
    rowN = img_input.shape[0]
    PR = int((1-np.sin(np.pi/3))*rowN/2)
    img_input = img_input[PR:rowN-PR, :]

    # set control parameters
    # closing operation kernel in module2
    Kernel1 = np.ones((3, 3))
    if filetype == 'MDI':
        # opening operation kernel in module2
        Kernel2 = np.ones((11, 11))
        # region gorwing threshold
        Thresh = 50
    elif filetype == 'HMI':
        Kernel2 = np.ones((9, 9))
        Thresh = 30
    # closing operation kernel in module4
    Kernel3 = np.ones((5, 5))
    # ARs area threshold in module4
    Size = 351
    # dilating operation kernel in module5
    Kernel4 = np.ones((23, 23))

    # detect ARs
    result = ARdetection(img_input, Kernel1, Kernel2,
                          Thresh, Kernel3, Size, Kernel4)

    # picture with ARs labeled
    img_label = result[6]
    num = np.max(img_label)

    CR_array = np.ones(num)*CR
    area = ARArea(img_label, img_input)
    flux = ARFlux(img_label, img_input)
    location = ARLocat(img_label, img_input)

    # parameters of all detected ARs
    param = np.zeros((12, num))
    param[0, :] = CR_array
    param[1, :] = np.linspace(1, num, num)
    param[2:8, :] = location
    param[8:10, :] = area
    param[10:12, :] = flux

    return param


# --------------------------------------------------------------------------------
if __name__ == '__main__':
    data_missing_file = []

    # batch process MDI and HMI magnetograms
    Filepath = 'D:\\python program\\活动区识别\\SynopMr'
    lists1 = os.listdir(Filepath)

    AR_parameters = np.zeros((12, 1))
    for j in range(len(lists1)):
        path = lists1[j]
        # folder of maps of each instrument 
        filepath_new = Filepath+'\\'+path
        lists2 = os.listdir(filepath_new)

        fileNum = len(lists2)
        print(fileNum)

        for i in range(fileNum):
            filename1 = lists2[i]
            # a single magnetogram
            file = os.path.join(filepath_new, filename1)

            # determine whether the data of map is missed
            if data_missing(file):
                # if data is missed, append the map to the list.
                data_missing_file.append(filename1)

            # use the length of file name to judge its source: HMI or MDI
            if len(filename1) > 35:
                filetype = 'HMI'
            else:
                filetype = 'MDI'

            parameters = GetARP(file, filetype)

            AR_parameters = np.concatenate((AR_parameters, parameters), axis=1)

            print(filename1+' done')
            print('      ')

        print(path.split()[1] + ' magnetogram done')

    # save the data
    AR_parameters = np.delete(AR_parameters, 0, axis=1)
    File2 = 'D:/python program/活动区识别/result data/Parameters of each AR(all).npy'
    np.save(File2, AR_parameters)
