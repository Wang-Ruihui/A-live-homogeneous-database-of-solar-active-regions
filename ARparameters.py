# -*- coding: utf-8 -*-
"""
Compute physical parameters of detected Active Regions (ARs)
"""

import numpy as np
from skimage import measure
from math import erf

# Solar constants
rs = 695.5  # Solar radius in Mm
microH = 2 * np.pi * rs**2 / 1e6  # One-millionth of a solar hemisphere area

# Synoptic magnetogram resolution
dx_deg_per_pixel = 360 / 3600  # Degrees per pixel in longitude
dy_sinlat_per_pixel = 2 / 1440  # Delta(sin(latitude)) per pixel in latitude

# Pixel area in physical units
pixel_area_Mm2 = (dx_deg_per_pixel / 180 * np.pi) * \
    dy_sinlat_per_pixel * rs**2  # Mm²
pixel_area_microH = pixel_area_Mm2 / microH  # In units of 1e-6 solar hemisphere


def GetCentroid(AR, locat):
    """
    Compute the flux-weighted centroid of an active region.

    Parameters
    ----------
    AR : array
        Magnetic field values (positive or negative) of pixels belonging to one polarity of an AR.
    locat : tuple of arrays
        (row_indices, col_indices) of the AR pixels in the full synoptic map.

    Returns
    -------
    lat : float
        Latitude of the centroid in degrees.
    lon : float
        Longitude of the centroid in degrees.
    """
    rowN, colN = 1440, 3600  # Image dimensions

    # Flux-weighted centroid in pixel coordinates
    row = np.sum(AR * locat[0]) / np.sum(AR)
    col = np.sum(AR * locat[1]) / np.sum(AR)

    # Convert to latitude and longitude
    sin_lat = (row + 0.5) / rowN * 2 - 1
    lat = np.arcsin(sin_lat) * 180 / np.pi
    lon = (col + 0.5) / colN * 360

    return lat, lon


def ARLocat(img_label, img_input):
    """
    Compute centroid positions for positive, negative, and total flux of each AR.

    Parameters
    ----------
    img_label : 2D array
        Labeled image where each AR has a unique integer label.
    img_input : 2D array
        Original synoptic magnetogram.

    Returns
    -------
    location : 2D array of shape (6, N_ARs)
        Each column corresponds to one AR:
        - Rows 0,1: latitude and longitude of positive polarity centroid
        - Rows 2,3: latitude and longitude of negative polarity centroid
        - Rows 4,5: latitude and longitude of total (unsigned) flux centroid
    """
    n_ars = int(np.max(img_label))
    location = np.zeros((6, n_ars))

    for i in range(1, n_ars + 1):
        mask = (img_label == i)
        AR = img_input[mask]
        pos_flux = AR[AR > 0]
        neg_flux = AR[AR < 0]

        # Total flux centroid
        location[4, i-1], location[5, i -
                                   1] = GetCentroid(np.abs(AR), np.where(mask))

        # Polarity-specific centroids
        pos_mask = (img_input > 0) & mask
        neg_mask = (img_input < 0) & mask

        if pos_flux.size > 0:
            location[0, i-1], location[1, i -
                                       1] = GetCentroid(pos_flux, np.where(pos_mask))
        if neg_flux.size > 0:
            location[2, i-1], location[3, i -
                                       1] = GetCentroid(np.abs(neg_flux), np.where(neg_mask))

    return location


def Distance(loc):
    """
    Compute great-circle angular distance between two points on a sphere.

    Parameters
    ----------
    loc : array-like, shape (4,)
        [lat0, lon0, lat1, lon1] in degrees.

    Returns
    -------
    delta : float
        Angular separation in degrees.
    """
    lat0, lon0, lat1, lon1 = np.deg2rad(loc)
    cos_delta = (
        np.cos(lat0) * np.cos(lat1) * np.cos(lon1 - lon0) +
        np.sin(lat0) * np.sin(lat1)
    )
    # Clamp to [-1, 1] to avoid numerical errors
    cos_delta = np.clip(cos_delta, -1.0, 1.0)
    delta = np.arccos(cos_delta)
    return np.rad2deg(delta)


def Tilt(loc):
    """
    Compute tilt angle of an active region bipole.

    Parameters
    ----------
    loc : array-like, shape (6,)
        [lat_pos, lon_pos, lat_neg, lon_neg, lat_center, lon_center] in degrees.

    Returns
    -------
    tilt : float
        Tilt angle in degrees (Joy's law convention).
    """
    lat0, lon0, lat1, lon1, latw, lonw = np.deg2rad(loc)
    dlat = lat0 - lat1
    dlon = (lon0 - lon1) * np.cos(latw)
    tilt = np.arctan2(-dlat, dlon)  # Use arctan2 for correct quadrant
    return np.rad2deg(tilt)


def ARArea(img_label, img_input):
    """
    Compute areas of positive and negative polarities for each AR.

    Returns
    -------
    area : array, shape (2, N_ARs)
        Row 0: positive polarity area (in micro-hemispheres)
        Row 1: negative polarity area
    """
    n_ars = int(np.max(img_label))
    area = np.zeros((2, n_ars))

    for i in range(n_ars):
        ar = img_input * (img_label == (i + 1))
        area[0, i] = np.sum(ar > 0) * pixel_area_microH  # Positive area
        area[1, i] = np.sum(ar < 0) * pixel_area_microH  # Negative area

    return area


def ARFlux(img_label, img_input):
    """
    Compute magnetic flux of each polarity for all ARs.

    Returns
    -------
    flux : array, shape (2, N_ARs)
        Magnetic flux in Maxwells (Mx).
        Row 0: positive flux, Row 1: negative flux.
    """
    n_ars = int(np.max(img_label))
    flux = np.zeros((2, n_ars))

    for i in range(1, n_ars + 1):
        AR = img_input[img_label == i]
        pos = AR[AR > 0]
        neg = AR[AR < 0]
        # 1 Gauss * Mm² = 1e16 Mx
        flux[0, i-1] = np.sum(pos) * pixel_area_Mm2 * 1e16
        flux[1, i-1] = np.sum(neg) * pixel_area_Mm2 * 1e16

    return flux


def ARBmax(img_label, img_input):
    """
    Compute maximum field strength within each AR.

    Returns
    -------
    Bmax : array, shape (N_ARs,)
        Peak |B| in Gauss.
    """
    n_ars = int(np.max(img_label))
    Bmax = np.zeros(n_ars)

    for i in range(1, n_ars + 1):
        AR = img_input[img_label == i]
        Bmax[i-1] = np.max(np.abs(AR))

    return Bmax


def FDF(img_label, img_input, lamr):
    """
    Compute the Final Dipole Field (FDF) contribution from each AR (Wang 2021).

    Parameters
    ----------
    lamr : float
        Supergranular diffusion scale (degrees).

    Returns
    -------
    fdf : array
        FDF values for all ARs.
    """
    n_ars, (M1, N1) = int(np.max(img_label)), img_label.shape
    dx = np.deg2rad(dx_deg_per_pixel)
    dy = dy_sinlat_per_pixel
    A = 0.21
    lamr_rad = np.deg2rad(lamr)

    fdf = []
    for i in range(n_ars):
        ar = np.zeros_like(img_input)
        mask = (img_label == i + 1)
        ar[mask] = img_input[mask]

        # Flux balancing (Yeates 2020)
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p, n = np.sum(pos), np.sum(np.abs(neg))
        if p > 0:
            ar[ar > 0] *= (p + n) / (2 * p)
        if n > 0:
            ar[ar < 0] *= (p + n) / (2 * n)

        fdf_val = 0.0
        for j in range(M1):
            lat = np.arcsin((j + 0.5) / M1 * 2 - 1)
            fdf_val += np.sum(ar[j, :]) * erf(lat / (lamr_rad * np.sqrt(2)))
        fdf_val = A * fdf_val * dx * dy
        fdf.append(fdf_val)

    return np.array(fdf)


def IDF(img_label, img_input):
    """
    Compute the Initial Dipole Field (IDF) contribution from each AR.

    Returns
    -------
    idf : array
        IDF values for all ARs.
    """
    n_ars, (M1, N1) = int(np.max(img_label)), img_label.shape
    dx = np.deg2rad(dx_deg_per_pixel)
    dy = dy_sinlat_per_pixel

    idf = []
    for i in range(n_ars):
        ar = np.zeros_like(img_input)
        mask = (img_label == i + 1)
        ar[mask] = img_input[mask]

        # Flux balancing (Yeates 2020)
        pos = ar[ar > 0]
        neg = ar[ar < 0]
        p, n = np.sum(pos), np.sum(np.abs(neg))
        if p > 0:
            ar[ar > 0] *= (p + n) / (2 * p)
        if n > 0:
            ar[ar < 0] *= (p + n) / (2 * n)

        idf_val = 0.0
        for j in range(M1):
            sinlat = (j + 0.5) / M1 * 2 - 1
            idf_val += np.sum(ar[j, :]) * sinlat
        idf_val = (3 / (4 * np.pi)) * idf_val * dy * dx
        idf.append(idf_val)

    return np.array(idf)


def FDFrs(flux, lat, sign, lamdr):
    """
    Compute FDF using empirical scaling laws (flux-latitude relation).

    Parameters
    ----------
    flux : array, shape (2, N_ARs)
        Magnetic flux of both polarities.
    lat : array, shape (N_ARs,)
        Latitude of AR center.
    sign : array, shape (N_ARs,)
        Sign of northern polarity (+1 or -1).
    lamdr : float
        Diffusion scale in degrees.

    Returns
    -------
    fdf_rs : array
        FDF values.
    """
    num = flux.shape[1]
    flux2 = np.max(np.abs(flux), axis=0)  # Use stronger polarity

    # Empirical separation-flux relation
    d_deg = (flux2 / 5.7e20) ** (1 / 1.73)  # Separation in degrees
    d_rad = np.deg2rad(d_deg)

    # Tilt angle (Joy's law)
    tilt_deg = np.abs(0.42 * lat)
    dlat_rad = d_rad * np.sin(np.deg2rad(tilt_deg))

    rs = 695.5  # Solar radius in Mm
    lat_rad = np.deg2rad(lat)
    idf_rs = sign * (3 / (4 * np.pi * rs**2)) * flux2 * \
        dlat_rad * np.cos(lat_rad) * 1e-16

    lamdr_rad = np.deg2rad(lamdr)
    a = np.sqrt(2 / np.pi) * (8 / 9)
    finf = a / lamdr_rad * np.exp(-lat_rad**2 / (2 * lamdr_rad**2))
    fdf_rs = finf * idf_rs

    return fdf_rs


def FDFBMR(flux, loc, lamdr):
    """
    Compute FDF using actual bipole geometry from detection.

    Parameters
    ----------
    flux : array, shape (2, N_ARs)
        Magnetic flux.
    loc : array, shape (6, N_ARs)
        AR locations from ARLocat().
    lamdr : float
        Base diffusion scale in degrees.

    Returns
    -------
    idf_bmr : array
        Initial dipole field.
    fdf_bmr : array
        Final dipole field.
    """
    lat = loc[4, :]
    lat_rad = np.deg2rad(lat)
    num = flux.shape[1]

    # Use mean unsigned flux (Yeates 2020)
    flux2 = np.mean(np.abs(flux), axis=0)

    # Latitude separation
    dlat_deg = loc[0, :] - loc[2, :]  # pos_lat - neg_lat
    dlat_rad = np.deg2rad(dlat_deg)

    rs = 695.5
    idf_bmr = (3 / (4 * np.pi * rs**2)) * flux2 * \
        dlat_rad * np.cos(lat_rad) * 1e-16

    # Effective diffusion scale
    d_sep_deg = Distance(loc[:4, :])
    lamdr_eff_deg = np.sqrt(lamdr**2 + (d_sep_deg / 2)**2)
    lamdr_eff_rad = np.deg2rad(lamdr_eff_deg)

    a = np.sqrt(2 / np.pi) * (8 / 9)
    finf = a / lamdr_eff_rad * np.exp(-lat_rad**2 / (2 * lamdr_eff_rad**2))
    fdf_bmr = finf * idf_bmr

    return idf_bmr, fdf_bmr
