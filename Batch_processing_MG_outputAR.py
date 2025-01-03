# -*- coding: utf-8 -*-
"""
batch process the magnetogram

生成活动区参数矩阵，保存到文件File2
"""

import numpy as np
from astropy.io import fits
import os
import shutil
from ARdetection import Get_ARi
from OutputARs import OutputARs
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


def check_and_clear_directory(directory_path):
    # check whether the path exist
    if not os.path.isdir(directory_path):
        print(f"路径 '{directory_path}' 不存在。")
        return

    # check whether the path is empty
    if not os.listdir(directory_path):
        print(f"目录 '{directory_path}' 已经为空。")
        return

    # clear the path
    for root, dirs, files in os.walk(directory_path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
    print(f"目录 '{directory_path}'清空完成。")

    return


# --------------------------------------------------------------------------------
if __name__ == '__main__':
    # output detected ARs

    data_missing_file = []
    # ARs in the same day
    ar_sameday = []
    # ARs in the map border
    ar_border = []

    # directory path of output ARs
    directory_path = 'D:\\test\\'
    # check and clear directory
    check_and_clear_directory(directory_path)

    # batch process the 23 and 24 solar cycle magneticgram
    Filepath = 'D:\\python program\\活动区识别\\SynopMr'
    # filepath of all instruments maps
    lists1 = os.listdir(Filepath)
    # map without ARs
    file_last = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'
    cr = 1909
    rt1 = 0.85
    rt2 = 1
    # the size of output AR map
    # size = (360, 180)
    size = (720, 360)

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

            # AR detection results
            results = Get_ARi(file, filetype, file_last,
                              file_next, rt1, rt2)

            img_org = results[0]
            img_label = results[8]

            # save AR maps
            # take the CR 1911 as initial contion, then CR 1912 is the initial time
            ar_sameday2, ar_border2 = OutputARs(
                directory_path, img_org, img_label, size, cr, 1912)
            if len(ar_sameday2) > 0:
                ar_sameday.append(ar_sameday2)
            if len(ar_border2) > 0:
                ar_border.append(ar_border2)

            cr += 1
            if cr == 1938:
                cr = 1941

            file_last = file

            print(filename1+' done')
            print('      ')

        print(path.split()[1] + ' magnetogram done')
