# -*- coding: utf-8 -*-
"""
Batch process magnetograms, output detected active region data, and save to files.
"""

import os
import numpy as np
import time
from ARdetection import Get_ARi
from Output_ARmap import output_ars
from Output_BoundaryARs import OB_main, merge_boundary_ARs


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


def process_magnetograms(output_dir, output_dir_boundaryAR, crs, cre, rt1=0.85, rt2=1,
                         size=(180, 360), lmax=300, type_b='KeepUF', allregion=False, RemoveRepeat=True):
    """
    Batch process magnetogram data and output active region maps.

    Parameters:
        output_dir (str): Output directory path.
        output_dir_boundaryAR (str): Output directory path for boundary ARs.
        rt1 (float): Parameter rt1, repeat AR threshold for decaying ARs
        rt2 (float): Parameter rt2, repeat AR threshold for emerging ARs
        size (tuple): Output image size.
        type_b (str): Method for balancing flux. Options are:
                      - 'KeepUF': Keep unbalanced flux.
                      - 'IncreSmall': Increase the smaller flux to match the larger one.
                      - 'DecreLarge': Decrease the larger flux to match the smaller one.

    Returns:
        date_ar(list): date of active regions.
        ar_border (list): Information about active regions at the boundary.
        ar_boundary (list): Information about active regions at the boundary.
    """
    date_ar = []
    ar_border = []
    ar_boundary = []

    for cr in range(crs, cre+1):
        # AR detection
        results = Get_ARi(cr, rt1=rt1, rt2=rt2, allregion=allregion,
                          RemoveRepeat=RemoveRepeat)
        img_org, img_label = results[0], results[8]

        # Save AR maps
        date_ar1, border = output_ars(
            output_dir, img_org, img_label, cr, 1912, size, lmax, type_b)
        date_ar.extend(date_ar1)
        ar_border.extend(border)

        ar_boundary1 = OB_main(cr, 1912, size, output_dir_boundaryAR)
        if ar_boundary1:
            ar_boundary.append(ar_boundary1)

        print(f"CR {cr} processed.")

    return date_ar, ar_border, ar_boundary


if __name__ == '__main__':
    start_time = time.time()

    output_dir_allAR = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\test2'
    output_dir_boundary = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\test2\boundary'

    crs = 2299  # 2225
    cre = 2303
    # Batch process magnetograms
    ar_sameday, ar_border, ar_boundary_all = process_magnetograms(
        output_dir_allAR, output_dir_boundary, crs, cre,
        type_b='KeepUF', allregion=False, RemoveRepeat=True)
    # allregion: True to detect all magentic region; False, regular AR detection

    existed_files, skipped_files = merge_boundary_ARs(
        output_dir_boundary, output_dir_allAR)

    print("All magnetograms processed.")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"程序运行时间: {elapsed_time:.4f} 秒")
