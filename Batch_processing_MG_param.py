# -*- coding: utf-8 -*-
"""
Batch process magnetograms to extract Active Region (AR) parameters.

This script reads synoptic magnetogram FITS files, detects active regions,
extracts their physical parameters, and saves them into a NumPy array and CSV file.
"""

import os
import re
import numpy as np
from astropy.io import fits
from ARdetection import Get_ARP


def has_missing_data(file_path):
    """
    Check if the magnetogram has missing data by inspecting NaN values near the equator.

    Parameters:
        file_path (str): Path to the FITS file.

    Returns:
        bool: True if missing data is detected, False otherwise.
    """
    try:
        with fits.open(file_path) as hdul:
            img_data = hdul[0].data
        # Check for NaN values in the equatorial region (row 720)
        nan_count = np.isnan(img_data[720, :]).sum()
        return nan_count > 0
    except Exception as e:
        print(f"[ERROR] Failed to read file {file_path}: {e}")
        return True


def get_cr_number_from_filename(filename):
    """
    Extract Carrington Rotation (CR) number from filename using regex.

    CR number is a 4-digit number in the filename.
    """
    match = re.search(r'(\d{4})', os.path.basename(filename))
    return int(match.group(1)) if match else None


def main():
    """
    Main function to batch process magnetogram files and extract AR parameters.
    """

    # Configuration
    output_file_npy = r'D:\python program\活动区识别\databseGithub\2025.10.8\database_allARs(lamr6.44)2.npy'
    output_file_csv = r'D:\python program\活动区识别\databseGithub\2025.10.8\database_allARs(lamr6.44)2.csv'

    lamr = 6.44     # Parameter used in AR detection algorithm

    crs = 1909
    cre = 2303

    # wheter remove repeat ARs
    RemoveRepeat = True
    print('RemoveRepeat is', RemoveRepeat)

    # List of files with missing data
    missing_data_files = []
    ar_param_list = []

    for cr in range(crs, cre+1):
        # Extract AR parameters
        try:
            params = Get_ARP(cr, lamr, allregion=False,
                             RemoveRepeat=RemoveRepeat)
            if params.size > 0:
                ar_param_list.append(params)
        except Exception as e:
            print(f"[ERROR] Error processing file {cr}: {e}")

        print(f"CR {cr} processed.")

    # Combine all AR parameter arrays horizontally
    if ar_param_list:
        ar_params_array = np.hstack(ar_param_list)

    else:
        print("[WARNING] No AR parameters extracted.")
        return

    # Save results
    np.save(output_file_npy, ar_params_array)
    np.savetxt(output_file_csv, ar_params_array.T, delimiter=',')
    print("Active Region parameters saved successfully.")

    # Output warning for files with missing data
    if missing_data_files:
        print("\n--- WARNING: Files with missing data ---")
        for f in missing_data_files:
            print(f"- {f}")


if __name__ == "__main__":
    main()
