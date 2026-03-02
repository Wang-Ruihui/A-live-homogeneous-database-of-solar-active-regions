# -*- coding: utf-8 -*-
"""
Created on Fri May 24 20:57:26 2024

output the detected ARs at the given size

@author: Ruihui Wang
"""
import numpy as np
import os
import time
from sunpy.coordinates.sun import carrington_rotation_time
from ARdetection import Get_ARi
from ARparameters import ARLocat
from Smooth_Resize_img import smooth_and_resize_img, SinlatToLat


def clear_directory(directory_path):
    """
    Clear all files in the specified directory.

    Parameters:
        directory_path (str): Directory path.

    Returns:
        None
    """
    if not os.path.isdir(directory_path):
        print(f"Path '{directory_path}' does not exist.")
        return

    for root, _, files in os.walk(directory_path, topdown=False):
        for file in files:
            os.remove(os.path.join(root, file))
    print(f"Directory '{directory_path}' cleared.")


def flux_balance(ar, type_b='KeepUF'):
    """
    Balance the flux of a 2D array (ar,180*360) in latitude and longitude grid.
    This method adjusts the positive and negative values of `ar` to balance their flux.

    Parameters:
        ar (np.ndarray): Input 2D array representing magnetic field values.
        type_b (str): Method for balancing flux. Options are:
                      - 'KeepUF': Keep unbalanced flux.
                      - 'IncreSmall': Increase the smaller flux to match the larger one.
                      - 'DecreLarge': Decrease the larger flux to match the smaller one.

    Returns:
        np.ndarray: The flux-balanced array.
    """
    # Ensure input is a NumPy array
    ar = np.asarray(ar).copy()
    m, n = ar.shape

    # Separate positive and negative components
    pos = np.where(ar > 0, ar, 0)  # Positive values
    neg = np.where(ar < 0, -ar, 0)  # Absolute values of negatives

    # Compute sine-weighted sums for each row
    sin_weights = np.sin((np.arange(m) + 0.5) * np.pi / m).reshape(-1, 1)
    p = np.sum(pos * sin_weights)  # Weighted sum of positive flux
    n = np.sum(neg * sin_weights)  # Weighted sum of negative flux

    # Apply the selected balancing method
    if type_b == 'KeepUF':
        # Scale both positive and negative flux to balance total flux
        scale_pos = (p + n) / (2 * p) if p > 0 else 1
        scale_neg = (p + n) / (2 * n) if n > 0 else 1
        ar[ar > 0] *= scale_pos
        ar[ar < 0] *= scale_neg

    elif type_b == 'IncreSmall':
        # Scale the smaller flux to match the larger one
        if p < n:
            ar[ar > 0] *= n / p if p > 0 else 1
        else:
            ar[ar < 0] *= p / n if n > 0 else 1

    elif type_b == 'DecreLarge':
        # Scale the larger flux to match the smaller one
        if p < n:
            ar[ar < 0] *= p / n if n > 0 else 1
        else:
            ar[ar > 0] *= n / p if p > 0 else 1

    else:
        raise ValueError(
            f"Invalid type_b: {type_b}. Choose from 'KeepUF', 'IncreSmall', or 'DecreLarge'.")

    return ar


def output_ars(directory_path, img_org, img_label, cr, cr0, size, lmax=300, type_b='KeepUF', arlat=None, latrange=(-90, 90)):
    """
    Output the magnetograms of detected ARs in the given size.

    Parameters:
        directory_path (str): Path to save the output files.
        img_org (np.ndarray): Original magnetograms.
        img_label (np.ndarray): Detected AR label of img.
        size (tuple): Size of the output image.
        cr (int): Carrington rotation number of the map.
        cr0 (int): Initial Carrington rotation number for time calculation.
        type_b (str): Method for balancing flux. Options are:
                      - 'KeepUF': Keep unbalanced flux.
                      - 'IncreSmall': Increase the smaller flux to match the larger one.
                      - 'DecreLarge': Decrease the larger flux to match the smaller one.
        arlat (list, optional): Latitudes of ARs. Defaults to None.
        latrange (list, optional): Latitude range of ARs to output. Defaults to [-90, 90].

    Returns:
        list: List of ARs appearing on the same day.
        list: List of ARs appearing on the border of the map.
    """
    nar = np.max(img_label)
    t0 = carrington_rotation_time(cr0)
    t1 = carrington_rotation_time(cr)
    t2 = carrington_rotation_time(cr + 1)

    # here give the resolution of original maps; (better if we keep the resolution unchanged during detection)
    row, col = img_org.shape

    date_counter = {}  # record the date of each AR emerging
    ar_border = []

    for i in range(nar):
        ar = np.zeros_like(img_org)
        index = (img_label == i + 1)
        ar[index] = img_org[index]

        # Calculate the date of the AR
        b1 = abs(ar[index])
        loc = np.where(index)
        col1 = np.sum(b1 * loc[1]) / np.sum(b1)
        lon = (col1 + 0.5)/col * 360

        # jduge the latitude of AR
        if arlat is not None:
            if not (latrange[0] <= arlat[i] <= latrange[1]):
                continue

        # remove the boundary ARs;they may be wrongly detected
        edge_lon_width = 20
        if lon > 360-edge_lon_width or lon < edge_lon_width:
            ar_border.extend([cr, i + 1])
            continue

        # Smooth and resize the AR map
        ar_lat = SinlatToLat(ar)
        ar_resized = smooth_and_resize_img(ar_lat, size, lmax)

        # Balance the flux
        ar_resized = flux_balance(ar_resized, type_b)

        # calculate the time of ar emerging
        t = t2 + lon / 360 * (t1 - t2)
        delta_day = round((t - t0).value)
        delta_day_str = f"{delta_day:05d}"

        # update the counter to get the number of AR emerging in the sameday
        if delta_day_str not in date_counter:
            date_counter[delta_day_str] = 0
        date_counter[delta_day_str] += 1
        n = date_counter[delta_day_str]  # 获取当前序号

        file_name = f'day_{delta_day_str}_{n}.txt'
        # Save the AR map
        file_path = os.path.join(directory_path, file_name)
        np.savetxt(file_path, ar_resized)  # 直接覆盖旧文件

    return date_counter, ar_border


if __name__ == '__main__':
    start_time = time.time()

    rt1 = 0.85
    rt2 = 1
    size = (180, 360)
    latrange = (-90, 0)  # Limit to the southern hemisphere
    directory_path = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\cycle25\1'

    date_counter = []
    ar_border = []

    crs = 2225
    cre = 2297

    for cr in range(crs, cre+1):
        results = Get_ARi(cr, rt1, rt2, RemoveRepeat=True)
        img_org, img_label = results[0], results[8]
        arlat = ARLocat(img_label, img_org)[4, :]

        # Save AR maps
        date_counter1, border = output_ars(
            directory_path, img_org, img_label, cr, 1912, size, lmax=300, type_b='KeepUF')
        date_counter.append(date_counter1)
        if border:
            ar_border.append(border)
        print(f"Processed CR: {cr}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"程序运行时间: {elapsed_time:.4f} 秒")
