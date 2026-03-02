# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:12:34 2024

Smooth the input image and resize to low resolution

@author: sky
"""

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import time


def SinlatToLat(imgSinlat):
    M, N = imgSinlat.shape
    eps_sinlat = 2/M*0.5
    sinlat_eq = np.linspace(-1+eps_sinlat, 1-eps_sinlat, M)
    eps_lat = np.pi/M*0.5
    latitudes = np.linspace(-np.pi/2+eps_lat, np.pi/2-eps_lat, M)

    # Interpolate each column independently
    f = interp1d(sinlat_eq, imgSinlat, axis=0, kind='linear',
                 bounds_error=False, fill_value="extrapolate")
    imgLat = f(np.sin(latitudes))

    return imgLat


def periodic_gaussian_filter1d(signal, sigma, period=360):
    """
    1D Gaussian filter with periodic boundary conditions.

    Parameters:
        signal (np.ndarray): Input 1D signal.
        sigma (float): Standard deviation of the Gaussian kernel.
        period (int): range to extend the boundary.

    Returns:
        np.ndarray: Smoothed signal with periodic boundaries handled.
    """
    # Extend the signal to handle circular convolution
    ext_signal = np.concatenate([signal[-period:], signal, signal[:period]])
    ext_smoothed = gaussian_filter1d(ext_signal, sigma, mode='nearest')
    return ext_smoothed[period:-period]


def adaptive_spherical_smoothing(data, lmax, print_sigma=False):
    """
    Apply latitude-adaptive Gaussian smoothing to spherical data.

    Parameters:
        data (np.ndarray): Input 2D array with shape (latitude, longitude)
        lmax (int): Maximum spherical harmonic degree. Controls the smoothing scale.
        print_sigma (bool): Whether to print the computed sigma values.

    Returns:
        np.ndarray: Smoothed data with reduced small-scale features.
    """
    smoothed = np.zeros_like(data)

    M_lat, N_lon = data.shape

    # Constant factor to convert FWHM to sigma in Gaussian filters
    cc = (8 * np.log(2)) ** 0.5

    # 1. Uniform smoothing along the latitude direction
    lat_res = 180 / M_lat  # Latitude resolution in degrees per pixel
    sigma_lat = 180 / (lmax * lat_res * cc)  # Convert to pixel-based sigma
    if print_sigma:
        print(f"Latitude smoothing parameter: σ={sigma_lat:.2f} pixels")

    for j in range(N_lon):
        smoothed[:, j] = gaussian_filter1d(
            data[:, j],
            sigma=sigma_lat,
            mode='nearest'  # Boundary handling
        )

    # 2. Longitude direction smoothing with latitude-dependent sigma
    latitudes = ((np.arange(M_lat) + 0.5) / M_lat * 180 - 90)
    lon_res_base = 360 / N_lon
    lon_res = np.cos(np.radians(latitudes)) * lon_res_base
    sigma_lon_lat = 180 / (lmax * lon_res * cc)

    # Pole protection
    sigma_lon_lat[np.abs(latitudes) > 85] = 180 / \
        (lmax * (np.cos(np.radians(85)) * lon_res_base) * cc)

    if print_sigma:
        print("Longitude smoothing done per row.")

    # Apply along each row using vectorization
    smoothed = np.vstack([
        periodic_gaussian_filter1d(row, sig, period=360)
        for row, sig in zip(smoothed, sigma_lon_lat)
    ])

    return smoothed


def smooth_and_resize_img(highres_data, size=(180, 360), lmax=300, print_sigma=False, eps=1e-10):
    """
    smooth and resized the high-resolution spherical magnetogram data to the low reolution.

    Steps:
        1. Latitude-adaptive Gaussian smoothing
        2. Spherical interpolation to a low-resolution grid

    Parameters:
        highres_data (np.ndarray): High-resolution input data (lat-lon grid)
        size (tuple (180,360)): Size of low-resolution data
        lmax (int): Maximum spherical harmonic degree used to control smoothing scale
        eps (float): Small value to avoid numerical issues at poles

    Returns:
        np.ndarray: Low-resolution output data after smoothing and interpolation
    """

    # 1. Apply adaptive spherical smoothing
    smoothed = adaptive_spherical_smoothing(highres_data, lmax, print_sigma)

    nlat_low, nlon_low = size
    nlat, nlon = highres_data.shape

    # Step 2: Define input and output grids
    lat_high = np.linspace(eps, np.pi - eps, nlat)
    lon_high = np.linspace(-np.pi + eps, np.pi - eps, nlon)

    lat_low = np.linspace(eps, np.pi - eps, nlat_low)
    lon_low = np.linspace(-np.pi + eps, np.pi - eps, nlon_low)

    # Create interpolator
    interpolating_function = RegularGridInterpolator(
        (lat_high, lon_high),
        smoothed,
        bounds_error=False,
        fill_value=0
    )

    # Create grid of low-res points
    lat_grid, lon_grid = np.meshgrid(lat_low, lon_low, indexing='ij')
    pts = np.column_stack((lat_grid.flatten(), lon_grid.flatten()))

    # Interpolate
    lowres_flat = interpolating_function(pts)
    lowres_data = lowres_flat.reshape(nlat_low, nlon_low)

    # Set small values to zero
    lowres_data[np.abs(lowres_data) < 1e-6] = 0

    return lowres_data


if __name__ == '__main__':
    start_time = time.time()

    # File path of input FITS file
    file_in = r'D:/python program/活动区识别/data/SynopticMap/HMI/HMI_polarfill/hmi.synoptic_mr_polfil_720s.2225.Mr_polfil.fits'

    # Load input FITS file
    input_file = fits.open(file_in)
    img_input = input_file[1].data
    input_file.close()

    # Convert sine-latitude projection to latitude-linear projection
    img_lat = SinlatToLat(img_input)

    # Set maximum spherical harmonic degree for smoothing
    lmax = 300
    # Perform preprocessing to get low-resolution magnetogram
    img_low_res = smooth_and_resize_img(
        img_lat, (180, 360), lmax, print_sigma=True)

    # Plot original and processed images
    fig, axes = plt.subplots(2, 1)
    axes[0].imshow(img_input, origin='lower', cmap='gray', vmin=-50, vmax=50)
    axes[1].imshow(img_low_res, origin='lower', cmap='gray', vmin=-50, vmax=50)
    plt.tight_layout()
    plt.show()

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"程序运行时间: {elapsed_time:.4f} 秒")
