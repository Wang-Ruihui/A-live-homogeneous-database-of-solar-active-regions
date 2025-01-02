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


# --------------------------------------------------------------------------------
if __name__ == '__main__':
    '''
    # 增加新磁图
    data_missing_file = []

    # batch process the 23 and 24 solar cycle magneticgram
    #Filepath = 'G:\\活动区识别\\SynopMr\\SynopMr HMI'
    Filepath = 'D:\\python program\\活动区识别\\SynopMr'
    lists1 = os.listdir(Filepath)

    AR_parameters = np.zeros((12, 1))

    # 2247-2097
    a = 2247-2097
    for i in range(a, len(lists1)):
        filename1 = lists1[i]
        file = os.path.join(Filepath, filename1)  # 每个综合磁图文件

        # 判断磁图低纬数据是否丢失
        if data_missing(file):
            # 如果数据缺失，标记，但继续识别
            data_missing_file.append(filename1)
            # continue
        continue

        # 判断文件类型，MDI 或 HMI
        if len(filename1) > 35:
            filetype = 'HMI'
        else:
            filetype = 'MDI'

        parameters = GetARP(file, filetype)

        AR_parameters = np.concatenate((AR_parameters, parameters), axis=1)
        # axis=1表示横向拼接

        print(filename1+' done')
        print('      ')

    AR_parameters = np.delete(AR_parameters, 0, axis=1)
    '''

    """
    #合并新磁图与旧磁图结果
    AR_parameters3 = np.concatenate((AR_parameters, AR_parameters2), axis=1)
    AR_parameters3 = np.transpose(AR_parameters3)
    np.savetxt('G:/活动区识别/data/allAR.csv', AR_parameters3, delimiter = ',')
    """

    # 全部MDI、HMI磁图识别
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
        # break

    # 保存数据
    AR_parameters = np.delete(AR_parameters, 0, axis=1)  # 删除第一列
    File2 = 'D:/python program/活动区识别/result data/Remove12/AR_Parameters_v1(lamr7.17,BL_yeates).npy'
    np.save(File2, AR_parameters)
