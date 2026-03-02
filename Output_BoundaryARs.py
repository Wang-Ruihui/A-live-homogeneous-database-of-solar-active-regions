# -*- coding: utf-8 -*-
"""
Process single magnetogram and detect Active Regions (ARs) in the boundary of map
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from sunpy.coordinates.sun import carrington_rotation_time

from ARdetection import ARdetection, ARdetection_org, ARdetection_allregion, Get_data, imshow
from ARparameters import ARLocat
from Output_ARmap import flux_balance
from Smooth_Resize_img import smooth_and_resize_img, SinlatToLat
import time


def get_kernels_and_threshold(cr):
    """
    Determine morphological kernels and thresholds based on Carrington Rotation number.
    HMI data (CR >= 2097) vs MDI data (CR < 2097).
    """
    if cr >= 2097:
        Kernel1 = np.ones((3, 3))   # CLOSE kernel
        Kernel2 = np.ones((9, 9))   # OPEN kernel
        Thresh = 30                 # region growing threshold
        Kernel4 = np.ones((19, 19))  # dilating kernel
        T_B = 100
    else:
        Kernel1 = np.ones((3, 3))
        Kernel2 = np.ones((11, 11))
        Thresh = 50
        Kernel4 = np.ones((19, 19))
        T_B = 130

    Kernel3 = np.ones((5, 5))       # closing kernel in module4
    Tarea = 351                      # area threshold for removing small regions

    return Kernel1, Kernel2, Thresh, Kernel3, Tarea, Kernel4, T_B


def load_images(cr):
    """
    Load synoptic maps for current CR and adjacent rotations (prev/next).
    Zero out overlapping longitudinal regions to avoid duplication.
    """
    img_current = Get_data(cr)
    img_prev = Get_data(cr - 1)
    img_prev2 = Get_data(cr - 2)
    img_next = Get_data(cr + 1)

    img_concat = np.concatenate((img_current, img_prev), axis=1)
    img_prev_concat = np.concatenate((img_prev, img_prev2), axis=1)
    img_next_concat = np.concatenate((img_next, img_current), axis=1)

    # Mask overlapping longitude zones
    for img in [img_concat, img_prev_concat, img_next_concat]:
        img[:, 0:3200] = 0
        img[:, 4000:7200] = 0

    return img_concat, img_prev_concat, img_next_concat


def output_ars_boundary(directory_path, file_name, img_org, size, lmax, type_b='KeepUF'):
    """
    Output detected AR magnetograms after processing.

    Parameters:
        directory_path (str): Output directory.
        file_name (str): Output filename.
        img_org (np.ndarray): Original AR magnetogram.
        size (tuple): Target output size (lat, lon).
        lmax (int): Maximum spherical harmonic degree (for smoothing).
        type_b (str): Flux balancing method:
                      - 'KeepUF': Keep unbalanced flux.
                      - 'IncreSmall': Increase smaller polarity to match larger.
                      - 'DecreLarge': Decrease larger polarity to match smaller.
    """
    ar = np.copy(img_org)
    ar_lat = SinlatToLat(ar)
    ar_resized = smooth_and_resize_img(ar_lat, size, lmax)
    ar_resized = flux_balance(ar_resized, type_b)

    file_path = os.path.join(directory_path, file_name)
    np.savetxt(file_path, ar_resized)


def extract_boundary_regions(output_dir, img_label, img_input, location, cr, cr0, size, lmax=300):
    """
    Extract ARs near the longitudinal boundary (0°/360°) and assign correct CR.
    Save them with filenames indicating emergence day relative to cr0.
    """
    lon_allars = location[5, :]
    edge_lon_width = 20
    indexes = np.where((lon_allars > 360 - edge_lon_width) &
                       (lon_allars < 360 + edge_lon_width))

    ar_boundary = []
    date_counter = {}

    for label_idx in indexes[0]:
        label = label_idx + 1
        mask = (img_label == label)
        ar = np.zeros_like(img_label)
        ar[mask] = img_input[mask]

        # Reconstruct boundary AR across 0° meridian
        ar_final = np.zeros_like(img_input[:, 0:3600])
        ar_final[:, 0:400] = ar[:, 3600:4000]
        ar_final[:, 3200:3600] = ar[:, 3200:3600]

        lon = location[5, label_idx]
        if lon > 360:
            cr_ar = cr - 1
            lon = lon - 360
        else:
            cr_ar = cr
        ar_boundary.append([cr_ar, label, lon])

        # Compute emergence day relative to reference time (cr0)
        t0 = carrington_rotation_time(cr0)
        t1 = carrington_rotation_time(cr_ar)
        t2 = carrington_rotation_time(cr_ar + 1)
        t = t2 + lon / 360 * (t1 - t2)
        delta_day = round((t - t0).value)
        delta_day_str = f"{delta_day:05d}"

        if delta_day_str not in date_counter:
            date_counter[delta_day_str] = 0
        date_counter[delta_day_str] += 1
        n = date_counter[delta_day_str]

        file_name = f'day_{delta_day_str}_{n}.txt'
        output_ars_boundary(output_dir, file_name, ar_final, size, lmax)

    return ar_boundary


def merge_boundary_ARs(input_dir, output_dir):
    """
    Merge all .txt boundary AR files from input_dir into output_dir.
    If a filename already exists, increment the suffix number (n) until unique.

    Parameters:
        input_dir (str): Input directory path.
        output_dir (str): Output directory path.

    Returns:
        tuple: (list of renamed files, list of skipped/failed files)
    """
    existed_files = []
    skipped_files = []

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if not filename.endswith('.txt'):
            continue

        src_path = os.path.join(input_dir, filename)
        base_name = filename[:10]  # e.g., 'day_09382_'

        try:
            arr = np.loadtxt(src_path)

            counter = 1
            new_filename = f"{base_name}{counter}.txt"
            dst_path = os.path.join(output_dir, new_filename)

            while os.path.exists(dst_path):
                counter += 1
                new_filename = f"{base_name}{counter}.txt"
                dst_path = os.path.join(output_dir, new_filename)
                existed_files.append(new_filename)

            np.savetxt(dst_path, arr)
            if new_filename != filename:
                print(f"Renamed file: {filename} -> {new_filename}")

        except Exception as e:
            print(f"[Error] Failed to process {filename}: {e}")
            skipped_files.append(filename)

    return existed_files, skipped_files


def plot_detection_result(output_path, img_input, img_border, location, num_ar):
    """
    Visualize AR detection results with contours and labels.

    Parameters:
        output_path (str): Path to save the figure (optional; currently shown only).
        img_input (np.ndarray): Input synoptic magnetogram.
        img_border (np.ndarray): Labeled AR boundary image for contouring.
        location (np.ndarray): AR location array from ARLocat().
        num_ar (int): Number of detected ARs.
    """
    Y, X = img_input.shape

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    ax.tick_params(length=4, width=1, labelsize=20)

    imshow(ax, img_input)
    ax.contour(img_border, colors='darkorange')

    ax.set_xticks(np.linspace(0, X, 9))
    ax.set_xticklabels(
        ['0', '90', '180', '270', '360', '90', '180', '270', '360'])
    ax.set_yticks(np.linspace(0, Y, 3))
    ax.set_yticklabels(
        ['sin(-90$^\circ$)', '0', 'sin(90$^\circ$)'], fontsize=16)

    ax.set_xlabel("Longitude (degree)", fontsize=24)

    lat_deg = location[4, :]
    lon_deg = location[5, :]

    lat_pixel = (np.sin(lat_deg / 180 * np.pi) +
                 np.sin(np.pi / 3)) / (np.sin(np.pi / 3) * 2) * Y
    lon_pixel = lon_deg * 10

    for i in range(num_ar):
        ax.text(lon_pixel[i] + 80, lat_pixel[i], i + 1, fontsize=18)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()


def OB_main(cr, cr0, size, output_dir_boundaryAR, allregion=False, RemoveRepeat=True):
    """
    Main pipeline for boundary AR detection in a given Carrington Rotation.

    Parameters:
        cr (int): Current Carrington Rotation number.
        cr0 (int): Reference CR for time indexing.
        size (tuple): Output AR map size.
        output_dir_boundaryAR (str): Directory to save boundary ARs.
        allregion (bool): If True, detect all regions without filtering.
        RemoveRepeat (bool): If True, remove duplicate ARs using temporal info.

    Returns:
        list or None: List of boundary ARs [cr_ar, label, lon], or None if failed.
    """
    if not os.path.exists(output_dir_boundaryAR):
        os.makedirs(output_dir_boundaryAR)

    try:
        # Load and preprocess images
        img0, img_last, img_next = load_images(cr)
        Kernel1, Kernel2, Thresh, Kernel3, Tarea, Kernel4, T_B = get_kernels_and_threshold(
            cr)

        # Detect ARs
        rt1, rt2 = 0.85, 1.0
        if allregion:
            _, _, _, _, result, img_label, img_border = ARdetection_allregion(
                img0, Kernel1, Kernel2, Thresh, Kernel3, Tarea, Kernel4)
        elif RemoveRepeat:
            _, _, _, _, _, _, result, img_label, img_border = ARdetection(
                img_last, img_next, rt1, rt2, img0, Kernel1, Kernel2,
                Thresh, Kernel3, Tarea, Kernel4, T_B)
        else:
            _, _, _, _, _, result, img_label, img_border = ARdetection_org(
                img0, Kernel1, Kernel2, Thresh, Kernel3, Tarea, Kernel4, T_B)

        num_ar = np.max(img_label)
        location = ARLocat(img_label, img0)

        # Optional: visualize results
        # plot_detection_result(None, img0, img_border, location, num_ar)

        # Extract and save boundary ARs
        ar_boundary = extract_boundary_regions(
            output_dir_boundaryAR, img_label, img0, location, cr, cr0, size, lmax=300)

        if ar_boundary:
            print(
                f'CR {cr}: processed, {len(ar_boundary)} boundary ARs detected.')
            return ar_boundary

    except Exception as e:
        print(f"[Error] Processing CR {cr}: {e}")
        return None


# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    start_time = time.time()

    output_dir_img = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\cycle25\test'
    output_dir_boundaryAR = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\cycle25\test\boundary'
    size = (180, 360)
    cr0 = 1912
    ar_boundary_all = []

    crs = 2299
    cre = 2303

    for cr in range(crs, cre + 1):
        ar_boundary1 = OB_main(
            cr, cr0, size, output_dir_boundaryAR, allregion=False, RemoveRepeat=True)
        if ar_boundary1:
            ar_boundary_all.append(ar_boundary1)
        print(f'{cr}: Processing completed')

    # Example usage of merge function (uncomment if needed):
    # output_dir_allAR = r'D:\WorkingProgram\sftprograms\ARmaps\Remove_repeat\cycle25\test'
    # existed_files, skipped_files = merge_boundary_ARs(output_dir_boundaryAR, output_dir_allAR)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total runtime: {elapsed_time:.4f} seconds")
