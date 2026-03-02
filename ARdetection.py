# -*- coding: utf-8 -*-
"""
Automatic detection of solar active regions from synoptic magnetograms

"""

import numpy as np
from collections import deque
import cv2
from pathlib import Path
from astropy.io import fits
from skimage import measure, morphology
from ARparameters import ARArea, ARFlux, ARLocat, ARBmax, IDF, FDF, FDFBMR

# ------------------------------------------------------------------------------


def imshow(ax, data, vmax=150, vmin=-150):
    """
    show image in the given range
    """
    ax.imshow(data, cmap='gray', vmin=vmin, vmax=vmax,
              origin='lower')  # cmap='gray'


def adaptiveThreshold(img_input, size, const):
    """
    adaptive threshold segmentation

    """
    img_smt = cv2.GaussianBlur(img_input, (size, size), 0, borderType=2)
    img_smt[img_smt < const] = const
    img_result = img_input - img_smt

    img_result[np.where(img_result > 0)] = 1
    img_result[np.where(img_result < 0)] = 0

    return img_result


def threshold(img_input, thresh=30):
    """
    intensity threshold segmentation
    """
    # thresh = 30
    img_thresh = (np.abs(img_input) > thresh).astype(np.uint8)

    # img_input1 = np.abs(img_input).astype("int16")
    # # img_input 'int16', param: 501,10
    # img_thresh = adaptiveThreshold(img_input1, 101, 10)

    return img_thresh

# ------------------------------------------------------------------------------


def region_grow(img, seeds, threshold, connectivity=8):
    """
    Region growing algorithm for image segmentation.

    Parameters:
        img (np.ndarray): Input 2D image array.
        seeds (list of tuples): List of seed points in (x, y) format.
        threshold (float or int): Threshold for similarity check.
        connectivity (int): 4 or 8, type of neighborhood connection.

    Returns:
        np.ndarray: Segmented binary mask with labeled regions.
    """
    h, w = img.shape
    seed_mask = np.zeros((h, w), dtype=np.uint8)
    label = 1  # 标记当前区域标签

    # 定义邻域偏移量
    if connectivity == 8:
        neighbors = [(-1, -1), (0, -1), (1, -1),
                     (1, 0), (1, 1), (0, 1),
                     (-1, 1), (-1, 0)]
    elif connectivity == 4:
        neighbors = [(0, -1), (1, 0), (0, 1), (-1, 0)]
    else:
        raise ValueError("connectivity must be 4 or 8")

    # 使用双端队列加速访问
    queue = deque()

    for x, y in seeds:
        if 0 <= x < h and 0 <= y < w and seed_mask[x, y] == 0:
            seed_mask[x, y] = label
            queue.append((x, y))

    while queue:
        x, y = queue.popleft()

        for dx, dy in neighbors:
            nx, ny = x + dx, y + dy
            if 0 <= nx < h and 0 <= ny < w:
                if seed_mask[nx, ny] == 0:
                    neighbor_val = abs(img[nx, ny])
                    if neighbor_val > threshold:
                        seed_mask[nx, ny] = label
                        queue.append((nx, ny))

    return seed_mask

# ------------------------------------------------------------------------------
# remove the repeat ARs


def LatRot(image, pixel_shifts):
    """
    Apply differential rotation to the image based on given pixel shifts for each latitude.

    Parameters:
    -----------
    image : np.ndarray
        Input image array (2D) representing solar magnetic field data.
    pixel_shifts : np.ndarray
        Array of integers specifying the number of pixels each row should be shifted. 
        Positive values shift right, negative left.

    Returns:
    --------
    np.ndarray
        The rotated image after applying differential rotation.
    """
    rows, cols = image.shape
    rotated_image = np.zeros((rows, cols))

    for lat_index in range(rows):
        shift = pixel_shifts[lat_index]

        if shift < 0:  # Shift to the left
            abs_shift = abs(shift)
            rotated_image[lat_index, 0:(cols - abs_shift)] \
                = image[lat_index, abs_shift:cols]
            rotated_image[lat_index, (cols - abs_shift):cols] \
                = image[lat_index, 0:abs_shift]
        elif shift > 0:  # Shift to the right
            rotated_image[lat_index,
                          shift:cols] = image[lat_index, 0:(cols - shift)]
            rotated_image[lat_index,
                          0:shift] = image[lat_index, (cols - shift):cols]
        else:  # No shift needed
            rotated_image[lat_index, :] = image[lat_index, :]

    return rotated_image


def RemoveRepeatAR(current_img, binary_img, prev_img, next_img, decayed_threshold, emerging_threshold):
    """
    Remove repeated Active Regions (ARs) by comparing with images from adjacent Carrington Rotations (CRs).

    Parameters:
    -----------
    current_img : np.ndarray
        Magnetic field image from the current CR.
    binary_img : np.ndarray
        Binary image after some processing step (module4).
    prev_img : np.ndarray
        Magnetic field image from the previous CR.
    next_img : np.ndarray
        Magnetic field image from the next CR.
    decayed_threshold : float
        Threshold for identifying repeat decayed ARs based on flux ratio.
    emerging_threshold : float
        Threshold for identifying not fully emerged ARs based on flux ratio.

    Returns:
    --------
    np.ndarray
        Processed binary image with repeated ARs removed.
    """
    rows, cols = current_img.shape
    processed_binary = np.copy(binary_img)

    # Mask out non-AR regions in the current image
    masked_current = current_img * binary_img

    # Label connected components in the binary image
    labeled_ar = measure.label(binary_img)

    # Differential rotation parameters
    omega_a = 13.562 - 13.20
    omega_b = -2.040
    omega_c = -1.487
    cr_duration_days = 27.27

    # Latitude sine values and corresponding angular velocities
    sin_latitudes = np.linspace(-1, 1, rows)
    angular_velocities = omega_a + omega_b * \
        sin_latitudes**2 + omega_c * sin_latitudes**4

    # Calculate pixel shifts for backward rotation (decay analysis)
    backward_pixel_shifts = (-angular_velocities *
                             cr_duration_days * 10).astype('int32')
    backward_rotated_img = LatRot(masked_current, backward_pixel_shifts)
    backward_rotated_labels = LatRot(
        labeled_ar, backward_pixel_shifts).astype('int64')

    # Analyze decayed ARs
    decayed_flux_ratios = np.zeros(np.max(backward_rotated_labels))
    decayed_polarities = np.zeros(np.max(backward_rotated_labels))
    analyze_ar(backward_rotated_img, prev_img, backward_rotated_labels,
               decayed_flux_ratios, decayed_polarities)

    # Calculate pixel shifts for forward rotation (emergence analysis)
    forward_pixel_shifts = (angular_velocities *
                            cr_duration_days * 10).astype('int32')
    forward_rotated_img = LatRot(masked_current, forward_pixel_shifts)
    forward_rotated_labels = LatRot(
        labeled_ar, forward_pixel_shifts).astype('int64')

    # Analyze emerging ARs
    emerging_flux_ratios = np.zeros(np.max(forward_rotated_labels))
    emerging_polarities = np.zeros(np.max(forward_rotated_labels))
    analyze_ar(forward_rotated_img, next_img, forward_rotated_labels,
               emerging_flux_ratios, emerging_polarities, dilate=True)

    # Identify and remove ARs that meet the criteria
    labels_to_remove = list(np.where((decayed_flux_ratios > decayed_threshold)
                                     & (decayed_polarities > -0.5))[0] + 1)
    labels_to_remove += list(np.where((emerging_flux_ratios > emerging_threshold)
                                      & (emerging_polarities > -0.5))[0] + 1)

    for label in labels_to_remove:
        mask = (labeled_ar == label)
        processed_binary[mask] = 0

    return processed_binary


def analyze_ar(rotated_img, comparison_img, rotated_labels, flux_ratios, polarities, dilate=False):
    """
    Helper function to analyze active regions based on flux ratios and polarities.

    Parameters:
    -----------
    rotated_img : np.ndarray
        Image after differential rotation.
    comparison_img : np.ndarray
        Image from another CR for comparison.
    rotated_labels : np.ndarray
        Labels of active regions after rotation.
    flux_ratios : np.ndarray
        Array to store calculated flux ratios.
    polarities : np.ndarray
        Array to store calculated polarities.
    dilate : bool, optional
        Whether to dilate the region for better matching due to evolution, default is False.
    """
    max_label = np.max(rotated_labels)
    for i in range(max_label):
        index = (rotated_labels == i + 1)
        if np.any(index):
            if dilate:
                index2 = cv2.dilate(index.astype('uint8'),
                                    np.ones((9, 9))).astype(bool)
            else:
                index2 = index
            new_region = rotated_img * index  # region after rotation
            old_region = comparison_img * index2  # region used for comparison

            new_abs_flux = np.sum(np.abs(new_region))
            new_net_flux = np.sum(new_region)
            old_abs_flux = np.sum(np.abs(old_region))
            old_net_flux = np.sum(old_region)

            if old_abs_flux > 0:
                flux_ratios[i] = old_abs_flux / new_abs_flux
                polarities[i] = np.sign(old_net_flux / new_net_flux) * abs(
                    new_net_flux / new_abs_flux) * abs(old_net_flux / old_abs_flux)


# ------------------------------------------------------------------------------

def RemoveDecayed(img_input, img_label, T_B, only_unipolar=False):
    """
    去除衰变活动区（包含单极区），并对剩余的活动区进行重新标记。

    参数:
    - img_input: 输入图像，包含活动区的信息。
    - img_label: 标签图像，每个活动区有唯一的标签值。
    - T_B: threshold for removing decayed ARs
    only_unipolar: True for only remove the unipolar regions; False for remove all decayed regions


    返回:
    - img_label: 更新后的标签图像，其中单极性活动区被移除。
    """
    max_label = np.max(img_label)
    flag = 1
    new_labels = np.zeros_like(img_label)  # 创建一个新的数组来存储新的标签

    for i in range(1, max_label + 1):
        mask = (img_label == i)
        AR = img_input[mask]

        pos = AR[AR > 0]
        neg = AR[AR < 0]

        if len(pos) == 0 or len(neg) == 0:
            continue  # 如果没有正或负磁通量，则跳过此区域

        meanB = np.mean(np.abs(AR))
        r_flux = np.sum(pos) / np.sum(np.abs(neg)
                                      ) if np.sum(np.abs(neg)) != 0 else float('inf')
        r_area = len(pos) / len(neg) if len(neg) != 0 else float('inf')

        con1 = (r_flux < 3) & (r_flux > 1/3)  # equal to |a-b|/(a+b) <0.5
        con2 = meanB > T_B
        # used to remove those regions with huge area imbalance
        con3 = (r_area < 3.5) & (r_area > 1/3.5)

        if only_unipolar:  # if only remove the unipolar region
            cons = con1
        else:
            cons = con1 & con2 & con3

        if cons:
            new_labels[mask] = flag
            flag += 1

    return new_labels


def ARdetection(img_last, img_next, rt1, rt2, img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4, T_B):
    """
    img_old  : image of the last CR
    rt1: threshold of flux ratio of the repeat decayed ARs
    rt2: threshold of flux ratio of the not-fully-emerged ARs
    img_input: input FITS file
    Kernel1  : closing operation kernel in module2
    Kernel2  : opening operation kernel in module2
    Thresh   : region grow threshold
    Kernel3  : closing operation kernel in module4
    Size     : the smallest size of AR,any region smaller than it will be removed.
    Kernel4  : dilating operation kernel in module5

    return:  result after each module,detection result, 
            boundary and labels of all ARs
    """
    img_input_copy = img_input.copy()
    plines = 130
    # remove the polar data
    img_input_copy[0:plines, :] = 0
    img_input_copy[1440-plines:1440, :] = 0

    # module1: intensity threshold segmentation, adaptive threshold
    module1 = threshold(img_input_copy, Thresh)

    # module2: morphological closing operation and opening operation to
    #          remove small magnetic segments and get the kernel pixels of ARs
    module2 = cv2.morphologyEx(module1, cv2.MORPH_CLOSE, Kernel1)
    module2 = cv2.morphologyEx(module2, cv2.MORPH_OPEN, Kernel2)

    # module3: region growing
    coords = np.column_stack(np.where(module2 > 0))
    seeds = [(x, y) for x, y in coords]
    module3 = region_grow(
        img_input, seeds, threshold=Thresh, connectivity=8)

    # module4: closing operation and remove small decayed ARs segments
    module4 = cv2.morphologyEx(module3, cv2.MORPH_CLOSE, Kernel3)
    module4 = morphology.remove_small_objects(module4.astype('bool'), Size)
    module4 = module4.astype("uint8")

    # # module42: removing the repeat ARs
    module42 = RemoveRepeatAR(img_input, module4, img_last, img_next, rt1, rt2)

    # module5: merging neighbor regions and removing unipolar regions
    # merge neighbor regions
    img_shape = cv2.dilate(module42, Kernel4)
    img_label = measure.label(img_shape)
    img_label = img_label * module42

    # remove unipolar regions
    img_label = RemoveDecayed(img_input, img_label, T_B)
    module5 = np.where(img_label >= 1, 1, 0)

    # get the border of each merged AR
    img_shape1 = module5.astype('uint8')
    img_shape2 = cv2.dilate(img_shape1, Kernel4)
    img_border = img_shape2 - cv2.erode(img_shape2, np.ones((3, 3)))

    # final detection result
    result = module5 * img_input

    return module1, module2, module3, module4, module42, module5, result, \
        img_label, img_border


def ARdetection_org(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4, T_B):
    """
    img_input: input FITS file
    Kernel1  : closing operation kernel in module2
    Kernel2  : opening operation kernel in module2
    Thresh   : region grow threshold
    Kernel3  : closing operation kernel in module4
    Size     : the smallest size of AR,any region smaller than it will be removed.
    Kernel4  : dilating operation kernel in module5

    return:  result after each module,detection result, 
            boundary and labels of all ARs
    """
    img_input_copy = img_input.copy()
    plines = 130
    # remove the polar data
    img_input_copy[0:plines, :] = 0
    img_input_copy[1440-plines:1440, :] = 0

    # module1: intensity threshold segmentation, adaptive threshold
    module1 = threshold(img_input_copy, thresh=Thresh)

    # module2: morphological closing operation and opening operation to
    #          remove small magnetic segments and get the kernel pixels of ARs
    module2 = cv2.morphologyEx(module1, cv2.MORPH_CLOSE, Kernel1)
    module2 = cv2.morphologyEx(module2, cv2.MORPH_OPEN, Kernel2)

    # module3: region growing
    coords = np.column_stack(np.where(module2 > 0))
    seeds = [(x, y) for x, y in coords]
    module3 = region_grow(
        img_input, seeds, threshold=Thresh, connectivity=8)

    # module4: closing operation and remove small decayed ARs segments
    module4 = cv2.morphologyEx(module3, cv2.MORPH_CLOSE, Kernel3)
    module4 = morphology.remove_small_objects(module4.astype('bool'), Size)
    module4 = module4.astype("uint8")

    # module5: merging neighbor regions and removing unipolar regions
    # merge neighbor regions
    img_shape = cv2.dilate(module4, Kernel4)
    img_label = measure.label(img_shape)
    img_label = img_label * module4

    # remove unipolar regions
    img_label = RemoveDecayed(img_input, img_label, T_B)
    module5 = np.where(img_label >= 1, 1, 0)

    # get the border of each merged AR
    img_shape1 = module5.astype('uint8')
    img_shape2 = cv2.dilate(img_shape1, Kernel4)
    img_border = img_shape2 - cv2.erode(img_shape2, np.ones((3, 3)))

    # final detection result
    result = module5 * img_input

    return module1, module2, module3, module4, module5, result, \
        img_label, img_border


def ARdetection_allregion(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4):
    """
    img_input: input FITS file
    Kernel1  : closing operation kernel in module2
    Kernel2  : opening operation kernel in module2
    Thresh   : region grow threshold
    Kernel3  : closing operation kernel in module4
    Size     : the smallest size of AR,any region smaller than it will be removed.
    Kernel4  : dilating operation kernel in module5

    return:  result after each module,detection result, 
            boundary and labels of all ARs
    """
    img_input_copy = img_input.copy()
    plines = 130
    # remove the polar data
    img_input_copy[0:plines, :] = 0
    img_input_copy[1440-plines:1440, :] = 0

    # module1: intensity threshold segmentation, adaptive threshold
    module1 = threshold(img_input_copy, Thresh)

    # module2: morphological closing operation and opening operation to
    #          remove small magnetic segments and get the kernel pixels of ARs
    module2 = cv2.morphologyEx(module1, cv2.MORPH_CLOSE, Kernel1)
    module2 = cv2.morphologyEx(module2, cv2.MORPH_OPEN, Kernel2)

    # module3: region growing
    coords = np.column_stack(np.where(module2 > 0))
    seeds = [(x, y) for x, y in coords]
    module3 = region_grow(
        img_input, seeds, threshold=Thresh, connectivity=8)

    # module4: closing operation and remove small decayed ARs segments
    module4 = cv2.morphologyEx(module3, cv2.MORPH_CLOSE, Kernel3)
    module4 = morphology.remove_small_objects(module4.astype('bool'), Size)
    module4 = module4.astype("uint8")

    img_label = measure.label(module4)
    # get the border of each merged AR
    img_border = cv2.dilate(module4, np.ones((7, 7))) - \
        cv2.dilate(module4, np.ones((5, 5)))

    # final detection result
    result = module4 * img_input

    return module1, module2, module3, module4, result, \
        img_label, img_border


# -------------------------------------------------------------------------------
'''
functions to detect ARs using the above detection modules
'''


def Get_data(cr, givenfile=None):
    """
    Read a FITS file of the input CR number and ensure the data is in float64 format.

    Parameters:
        cr (int): the CR number of the FITS file.

    Returns:
        np.ndarray: Data from the FITS file in float64 format.
    """
    if givenfile:
        file_path = givenfile
    else:
        if cr < 2097:
            file_path = f'D:/python program/活动区识别/SynopMr/SynopMr 1MDI/mdi.synoptic_Mr_96m.{cr}.data.fits'
        else:
            file_path = f'D:/python program/活动区识别/SynopMr/SynopMr HMI/hmi.Synoptic_Mr_720s.{cr}.synopMr.fits'

    if not Path(file_path).exists():
        print(f"File not found: {file_path},retun zero array")
        return np.zeros((1440, 3600))

    # Check if there are NaN values near the equator of the magnetogram (indicating missing data).
    with fits.open(file_path) as hdul:
        data = hdul[0].data  # 读取主数据
        equator_row = data.shape[0] // 2
        if np.isnan(data[equator_row, :]).any():
            print(f"There are missing data in the file: {file_path}.")
    img_input = data.astype(np.float64)

    # remove the nan value
    img_input[np.isnan(img_input)] = 0

    return img_input


def Get_ARi(cr, rt1=0.85, rt2=1, givenfile=None, file_type=None, allregion=False, RemoveRepeat=True):
    """
    # get AR detection image of File after each module

    cr: the Carrington rotation number of image
    rt1: threshold of flux ratio of the repeat decayed ARs
    rt2: threshold of flux ratio of the not-fully-emerged ARs
    RemoveRepeat: True for removing the repeat ARs; False for keeping

    givenfile and file_type are used to process the given files by users
    """

    img_input = Get_data(cr, givenfile)

    # set control parameters
    # closing operation kernel in module2
    Kernel1 = np.ones((3, 3))

    if not file_type:
        if cr < 2097:
            file_type = 'MDI'
        else:
            file_type = 'HMI'

    if file_type == 'MDI':
        # MDI synoptic magnetogram
        # opening operation kernel in module2
        Kernel2 = np.ones((11, 11))
        # region gorwing threshold
        Thresh = 50
        # dilating operation kernel in module5
        Kernel4 = np.ones((19, 19))
        # threshold to remove weak decayed ARs
        T_B = 130

    elif file_type == 'HMI':
        # HMI synoptic magnetogram
        Kernel2 = np.ones((9, 9))  # OPEN kernel
        Thresh = 30  # region gorwing threshold
        # dilating operation kernel in module5
        Kernel4 = np.ones((19, 19))
        # threshold to remove weak decayed ARs
        T_B = 100
    else:
        print('file_type should be MDI or HMI')

    # closing operation kernel in module4
    Kernel3 = np.ones((5, 5))  # Opening and closing kernel
    # ARs area threshold in module4
    Size = 351

    if allregion:
        module1, module2, module3, module4, result, img_label, img_border\
            = ARdetection_allregion(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4)
        module42 = 0
        module5 = 0
    elif RemoveRepeat:
        img_last = Get_data(cr-1)  # imgs used to remove repeat ARs
        img_next = Get_data(cr+1)
        # detection with removing repeat ARs
        module1, module2, module3, module4, module42, module5, result, img_label, img_border\
            = ARdetection(img_last, img_next, rt1, rt2, img_input,
                          Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4, T_B)
    else:
        # origional detection, keep repeat ARs
        module1, module2, module3, module4, module5, result, img_label, img_border\
            = ARdetection_org(img_input, Kernel1, Kernel2, Thresh, Kernel3, Size, Kernel4, T_B)
        module42 = 0

    return img_input, module1, module2, module3, module4, module42, module5, \
        result, img_label, img_border


def Get_ARP(cr, lamr, rt1=0.85, rt2=1, givenfile=None, file_type=None, allregion=False, RemoveRepeat=True):
    """
    get AR parameters from synoptic magnetograms

    cr: the Carrington rotation number of image
    rt1: threshold of flux ratio of the repeat decayed ARs
    rt2: threshold of flux ratio of the not-fully-emerged ARs
    lamr: value of the dynamo effectivity range

    Returns：param (17*num)(num: AR number in this CR)

           param[0,:]: Carrington rotation
           param[2:4,:]: position of the AR positive polarity
           param[4:6,:]: position of the AR negative polarity
           param[6:8,:]: position of the whole AR
           param[8:10,:]: area of the AR positive polarity and negative polarity
           param[10,:]: flux of the AR positive polarity
           param[11,:]: flux of the AR negative polarity   
           param[12,:]: max unsigned magnetic field of AR             
           param[13, :] = idf
           param[14, :] = fdf    
           param[15, :] = idf with BMR approximation
           param[16, :] = fdf with BMR approximation
           for all position, the first row is latitude and the second is longitude
           parameters in each column is for an same AR
    """
    # extraction ARs from magnet0gram
    result = Get_ARi(cr, rt1, rt2, givenfile,
                     file_type, allregion, RemoveRepeat)

    # segment picture,result[4] is the module5,the result of Closing operation
    img_input = result[0]
    img_label = result[8]
    num = np.max(img_label)

    CR_array = np.ones(num)*cr
    area = ARArea(img_label, img_input)
    flux = ARFlux(img_label, img_input)
    location = ARLocat(img_label, img_input)
    # the max unsigned magnetic field of AR
    Bmax = ARBmax(img_label, img_input)

    idf = IDF(img_label, img_input)
    fdf = FDF(img_label, img_input, lamr)
    # the BMR approximation of idf and fdf
    idfbmr, fdfbmr = FDFBMR(flux, location, lamr)
    # fdfsm = FDFsm(img_label, img_input, lamr)

    param = np.zeros((17, num))
    param[0, :] = CR_array
    param[1, :] = np.linspace(1, num, num)
    param[2:8, :] = location
    param[8:10, :] = area
    param[10:12, :] = flux
    param[12, :] = Bmax
    param[13, :] = idf
    param[14, :] = fdf
    param[15, :] = idfbmr
    param[16, :] = fdfbmr

    return param
