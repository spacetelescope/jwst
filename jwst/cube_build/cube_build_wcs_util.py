"""Routines related to WCS procedures of cube_build."""

import numpy as np
from jwst.assign_wcs import nirspec
from jwst.assign_wcs.util import wrap_ra
from gwcs import wcstools
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# ******************************************************************************
def find_corners_miri(input_data, this_channel, instrument_info, coord_system):
    """
    For MIRI channel data find the footprint of this data on the sky.

    For a specific channel on an exposure find the min and max of the
    spatial coordinates, either in alpha,beta or ra,dec depending
    on the type of cube being build. Also find the min and max of
    wavelength this channel covers.

    Parameters
    ----------
    input_data : IFUImage model
       Input model (or file)
    this_channel : str
       Channel working with
    instrument_info : dict
       Dictionary holding x pixel min and max values for each channel
    coord_system : str
       Coordinate system of output cube: skyalign, ifualign, internal_cal

    Returns
    -------
    a_min : float
        Minimum value of coord 1 - along axis 1 of the IFUcube
    b1 : float
        Coord 2 value corresponding  to a_min
    a_max : float
        Maximum value of coord 1 - along the axis 1 of the IFUcube
    b2 : float
        Coord 2 value corresponding to a_min
    a1 : float
        Coord 1 value coorsponding to b_min
    b_min : float
        Minimum value of coord 2 - along the axis 2 of the IFU cube
    a2 : float
        Coord 1 value coorsponding to b_max
    b_max : float
        Maximum value of coord 2 - along the axis 2 of the IFU cube
    lambda_min : float
        Minimum wavelength
    lambda_max : float
        Maximum wavelength
    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength

    xstart, xend = instrument_info.get_miri_slice_endpts(this_channel)
    ysize = input_data.data.shape[0]

    y, x = np.mgrid[:ysize, xstart:xend]

    if coord_system == "internal_cal":
        # coord1 = along slice
        # coord2 = across slice
        detector2alpha_beta = input_data.meta.wcs.get_transform("detector", "alpha_beta")
        coord1, coord2, lam = detector2alpha_beta(x, y)

        valid = np.logical_and(np.isfinite(coord1), np.isfinite(coord2))
        coord1 = coord1[valid]
        coord2 = coord2[valid]
        coord1 = coord1.flatten()
        coord2 = coord2.flatten()
        lam = lam.flatten()

        xedge = x + 1
        yedge = y + 1
        valid = np.where(yedge < 1023)
        xedge = xedge[valid]
        yedge = yedge[valid]
        # lam ~0 for this transform
        coord2b, coord1b, lamb = detector2alpha_beta(xedge, yedge)
        valid = np.logical_and(np.isfinite(coord1b), np.isfinite(coord2b))
        coord1b = coord1b[valid]
        coord2b = coord2b[valid]
        coord1b = coord1b.flatten()
        coord2b = coord2b.flatten()
        lamb = lamb.flatten()
        coord1 = np.concatenate([coord1, coord1b])
        coord2 = np.concatenate([coord2, coord2b])
        lam = np.concatenate([lam, lamb])

    else:  # skyalign or ifualign
        # coord1 = ra
        # coord2 = dec
        coord1, coord2, lam = input_data.meta.wcs(x, y)
        valid = np.logical_and(np.isfinite(coord1), np.isfinite(coord2))
        coord1 = coord1[valid]
        coord2 = coord2[valid]
        # fix 0/360 wrapping in ra. Wrapping makes it difficult to determine
        # ra range
        coord1_wrap = wrap_ra(coord1)
        coord1 = coord1_wrap

    coord1 = coord1.flatten()
    coord2 = coord2.flatten()
    a_min = np.nanmin(coord1)
    a_max = np.nanmax(coord1)

    # find index of min a value
    a1_index = np.argmin(coord1)
    a2_index = np.argmax(coord1)

    b1 = coord2[a1_index]
    b2 = coord2[a2_index]

    b_min = np.nanmin(coord2)
    b_max = np.nanmax(coord2)

    # find index of min b value
    b1_index = np.argmin(coord2)
    b2_index = np.argmax(coord2)
    a1 = coord1[b1_index]
    a2 = coord1[b2_index]

    lambda_min = np.nanmin(lam)
    lambda_max = np.nanmax(lam)

    if coord_system != "internal_cal":
        # before returning,  ra should be between 0 to 360
        a_min %= 360
        a_max %= 360

        a1 %= 360
        a2 %= 360

    return a_min, b1, a_max, b2, a1, b_min, a2, b_max, lambda_min, lambda_max


# *****************************************************************************


def find_corners_nirspec(input_data, coord_system):
    """
    Find the sky footprint of a slice of a NIRSpec exposure.

    For each slice find:
    a. the min and max spatial coordinates (along slice, across slice) or
       (ra,dec) depending on coordinate system of the output cube.
    b. min and max wavelength

    Parameters
    ----------
    input_data : IFUImageModel
       Input calibrated model (or file)
    coord_system : str
       Coordinate system of output cube: skyalign, ifualign, internal_cal

    Returns
    -------
    a_min : float
        Minimum value of coord 1 - along axis 1 of the IFUcube
    b1 : float
        Coord 2 value corresponding  to a_min
    a_max : float
        Maximum value of coord 1 - along the axis 1 of the IFUcube
    b2 : float
        Coord 2 value corresponding to a_min
    a1 : float
        Coord 1 value coorsponding to b_min
    b_min : float
        Minimum value of coord 2 - along the axis 2 of the IFU cube
    a2 : float
        Coord 1 value coorsponding to b_max
    b_max : float
        Maximum value of coord 2 - along the axis 2 of the IFU cube
    lambda_min : float
        Minimum wavelength
    lambda_max : float
        Maximum wavelength
    """
    nslices = 30
    a_slice = np.zeros(nslices * 2)
    b = np.zeros(nslices * 2)
    b_slice = np.zeros(nslices * 2)
    a = np.zeros(nslices * 2)
    lambda_slice = np.zeros(nslices * 2)
    k = 0
    # for NIRSPEC there are 30 regions
    log.info("Looping over slices to determine cube size")

    for i in range(nslices):
        slice_wcs = nirspec.nrs_wcs_set_input(input_data, i)
        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
        if coord_system == "internal_cal":
            # coord1 = along slice
            # coord2 = across slice
            detector2slicer = slice_wcs.get_transform("detector", "slicer")
            coord2, coord1, lam = detector2slicer(x, y)  # lam ~0 for this transform
            valid = np.logical_and(np.isfinite(coord1), np.isfinite(coord2))
            coord1 = coord1[valid]
            coord2 = coord2[valid]

            lmax = np.nanmax(lam.flatten())
            if lmax < 0.0001:
                lam = lam[valid] * 1.0e6
            coord1 = coord1.flatten()
            coord2 = coord2.flatten()
            lam = lam.flatten()
        else:  # coord_system: skyalign, ifualign
            # coord1 = ra
            # coord2 = dec
            coord1, coord2, lam = slice_wcs(x, y)
            valid = np.logical_and(np.isfinite(coord1), np.isfinite(coord2))
            coord1 = coord1[valid]
            coord2 = coord2[valid]
            lam = lam[valid]
            # fix 0/360 wrapping in ra. Wrapping makes it difficult to determine
            # ra range
            coord1_wrap = wrap_ra(coord1)
            coord1 = coord1_wrap
        # ________________________________________________________________________________
        coord1 = coord1.flatten()
        coord2 = coord2.flatten()
        lam = lam.flatten()

        a_slice[k] = np.nanmin(coord1)
        a_slice[k + 1] = np.nanmax(coord1)
        a1_index = np.argmin(coord1)
        a2_index = np.argmax(coord1)
        b1 = coord2[a1_index]
        b2 = coord2[a2_index]
        b[k] = b1
        b[k + 1] = b2

        b_slice[k] = np.nanmin(coord2)
        b_slice[k + 1] = np.nanmax(coord2)

        b1_index = np.argmin(coord2)
        b2_index = np.argmax(coord2)
        a1 = coord1[b1_index]
        a2 = coord1[b2_index]
        a[k] = a1
        a[k + 1] = a2

        lambda_slice[k] = np.nanmin(lam)
        lambda_slice[k + 1] = np.nanmax(lam)

        k = k + 2
    # ________________________________________________________________________________
    # now test the ra slices for consistency. Adjust if needed.

    a_min = np.nanmin(a_slice)
    a_max = np.nanmax(a_slice)
    a1_index = np.argmin(a_slice)
    a2_index = np.argmax(a_slice)
    b1 = b[a1_index]
    b2 = b[a2_index]

    b_min = np.nanmin(b_slice)
    b_max = np.nanmax(b_slice)
    b1_index = np.argmin(b_slice)
    b2_index = np.argmax(b_slice)
    a1 = a[b1_index]
    a2 = a[b2_index]

    lambda_min = min(lambda_slice)
    lambda_max = max(lambda_slice)

    return a_min, b1, a_max, b2, a1, b_min, a2, b_max, lambda_min, lambda_max
