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
    用于判断文件中数据是否丢失，如果磁图赤道附近有nan值，则有数据丢失

    Parameters
    ----------
    file : 文件地址，fits文件

    Returns
    -------
    bool
        文件数据丢失返回True，没有丢失返回False
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
    # 检查目录是否存在
    if not os.path.isdir(directory_path):
        print(f"路径 '{directory_path}' 不存在。")
        return

    # 检测目录是否为空
    if not os.listdir(directory_path):
        print(f"目录 '{directory_path}' 已经为空。")
        return

    # 清空目录
    for root, dirs, files in os.walk(directory_path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
    print(f"目录 '{directory_path}'清空完成。")

    return


# --------------------------------------------------------------------------------
if __name__ == '__main__':
    '''
    # 增加新磁图
    data_missing_file = []

    # batch process the 23 and 24 solar cycle magneticgram
    #Filepath = 'G:\\活动区识别\\SynopMr\\SynopMr HMI'
    Filepath = 'D:\\python program\\活动区识别\\SynopMr\\SynopMr HMI'
    lists1 = os.listdir(Filepath)

    # 2247-2097
    cr = 2266
    rt = 0.7
    a = cr-2097
    file_old = 'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.' + \
        str(cr-1)+'.synopMr.fits'
    for i in range(a, len(lists1)):
        filename1 = lists1[i]
        file = os.path.join(Filepath, filename1)  # 每个综合磁图文件

        # 判断磁图低纬数据是否丢失
        if data_missing(file):
            # 如果数据缺失，标记，但继续识别
            data_missing_file.append(filename1)
            # continue

        # 判断文件类型，MDI 或 HMI
        if len(filename1) > 35:
            filetype = 'HMI'
        else:
            filetype = 'MDI'

        # AR detection results
        results = Get_ARi(file, filetype, file_old, rt)
        # transform to uint8 to reduce the size
        img_label = results[8].astype("uint8")

        # save the img_label
        File2 = 'D:/python program/活动区识别/result data/img label(5,rt0.7,close)/CR ' + \
            str(cr)+'.npy'
        np.save(File2, img_label)
        cr += 1

        file_old = file

        print(filename1+' done')
        print('      ')


    '''
    # output detected ARs
    # 全部MDI、HMI磁图识别
    data_missing_file = []
    ar_sameday = []
    ar_border = []

    # directory path of output ARs
    directory_path = 'D:\\WorkingProgram\\sftprograms\\ARmaps\\Remove_NotEmerge_Repeat\\dilate19\\balancedY(360,720)\\'
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
    # size = (360, 180)
    size = (720, 360)

    for j in range(len(lists1)):
        path = lists1[j]
        filepath_new = Filepath+'\\'+path  # 每个仪器的综合磁图的文件夹
        lists2 = os.listdir(filepath_new)

        fileNum = len(lists2)
        print(fileNum)

        for i in range(fileNum):
            filename1 = lists2[i]
            file = os.path.join(filepath_new, filename1)  # 每个综合磁图文件
            if i+1 < fileNum:
                file_next = os.path.join(filepath_new, lists2[i+1])
            elif j == 0:
                file_next = 'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.2097.synopMr.fits'
            else:
                file_next = 'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.2077.data.fits'

            # 判断磁图低纬数据是否丢失
            if data_missing(file):
                # 如果数据缺失，标记，但继续识别
                data_missing_file.append(filename1)
                # continue

            # 判断文件类型，MDI 或 HMI
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
            # take the CR1911 as initial contion, then C1912 is the initial time
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
