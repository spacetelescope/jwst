""" Routines related to WCS procedures of cube_build
"""
import numpy as np
from ..assign_wcs import nirspec
from gwcs import wcstools
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# ********************************************************************************
def find_corners_MIRI(input, this_channel, instrument_info, coord_system):

    """ For MIRI channel data find the foot of this data on the sky

    For a specific channel on an exposure find the min and max of the
    spatial coordinates, either in  alpha,beta or ra,dec dedpending
    on the type of cube being build. Also find the min and max of
    wavelength this channel covers.

    Parameters
    ----------
    input : data model
       input model (or file)
    this_channel : str
       channel working with
    instrument_info : dictionary
       dictionary holding x pixel min and max values for each channel
    coord_system : str
       coordinate system of output cube: skyalign, ifualign, internal_cal

    Returns
    -------
    corners of the spatial system determined from min and max of spatial values
    min and max of wavelength

    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength

    xstart, xend = instrument_info.GetMIRISliceEndPts(this_channel)
    ysize = input.data.shape[0]

    y, x = np.mgrid[:ysize, xstart:xend]

    if coord_system == 'internal_cal':
        # coord1 = along slice
        # coord2 = across slice
        detector2alpha_beta = input.meta.wcs.get_transform('detector',
                                                           'alpha_beta')
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
        coord2b, coord1b, lamb = detector2alpha_beta(xedge,yedge) #lam ~0 for this transform
        valid = np.logical_and(np.isfinite(coord1b), np.isfinite(coord2b))
        coord1b = coord1b[valid]
        coord2b = coord2b[valid]
        coord1b = coord1b.flatten()
        coord2b = coord2b.flatten()
        lamb = lamb.flatten()
        coord1 = np.concatenate([coord1,coord1b])
        coord2 = np.concatenate([coord2,coord2b])
        lam = np.concatenate([lam,lamb])

    else:  # skyalign or ifualign
        # coord1 = ra
        # coord2 = dec
        coord1, coord2, lam = input.meta.wcs(x, y)
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

    a1_index = np.argmin(coord1)
    a2_index = np.argmax(coord1)

    b1 = coord2[a1_index]
    b2 = coord2[a2_index]

    b_min = np.nanmin(coord2)
    b_max = np.nanmax(coord2)

    b1_index = np.argmin(coord2)
    b2_index = np.argmax(coord2)
    a1 = coord1[b1_index]
    a2 = coord1[b2_index]

    lambda_min = np.nanmin(lam)
    lambda_max = np.nanmax(lam)

    return a_min, b1, a_max, b2, a1, b_min,  a2, b_max, lambda_min, lambda_max
# ********************************************************************************

def find_corners_NIRSPEC(input, this_channel, instrument_info, coord_system):

    """For a NIRSPEC slice on an exposure find the footprint of this data on the sky

    For each slice find:
    a. the min and max spatial coordinates (along slice, across slice) or (ra,dec) depending
       on coordinate system of the output cube.
    b. min and max wavelength

    Parameters
    ----------
    input: data model
       input model (or file)
    coord_system : str
       coordinate system of output cube: skyalign, ifualign, internal_cal

    Notes
    -----
    Returns
    -------
    min and max spaxial coordinates and wavelength for slice.

    """
    nslices = 30
    a_slice = np.zeros(nslices * 2)
    b = np.zeros(nslices * 2)
    b_slice = np.zeros(nslices * 2)
    a = np.zeros(nslices * 2)
    lambda_slice = np.zeros(nslices * 2)
    k = 0
    # for NIRSPEC there are 30 regions
    log.info('Looping over slices to determine cube size .. this takes a while')

    for i in range(nslices):
        slice_wcs = nirspec.nrs_wcs_set_input(input, i)
        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
        if coord_system == 'internal_cal':
            # coord1 = along slice
            # coord2 = across slice
            detector2slicer = slice_wcs.get_transform('detector','slicer')
            coord2, coord1, lam = detector2slicer(x,y) #lam ~0 for this transform
            valid = np.logical_and(np.isfinite(coord1), np.isfinite(coord2))
            coord1 = coord1[valid]
            coord2 = coord2[valid]

            #TODO - fix kludge after figure out how to determine what the units
            # for wavelength are coming out of detector2slicer.
            #print(input.meta.wcs.output_frame.unit[2])
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
        a_slice[k + 1] =np.nanmax(coord1)
        a1_index = np.argmin(coord1)
        a2_index = np.argmax(coord1)
        b1 = coord2[a1_index]
        b2 = coord2[a2_index]
        b[k] = b1
        b[k+1] = b2

        b_slice[k] = np.nanmin(coord2)
        b_slice[k + 1] = np.nanmax(coord2)

        b1_index = np.argmin(coord2)
        b2_index = np.argmax(coord2)
        a1 = coord1[b1_index]
        a2 = coord1[b2_index]
        a[k] = a1
        a[k +1] = a2

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
# _______________________________________________________________________________


def wrap_ra(ravalues):
    """Test for 0/360 wrapping in ra values.

    If exists it makes it difficult to determine
    ra range of IFU cube. So put them all on "one side" of 0/360 border

    Parameters
    ----------
    input : ravalues
        RA values numpy.ndarray

    Returns
    ------
    a numpy array of ra values all on "same side" of 0/360 border
    """

    valid = np.isfinite(ravalues)
    index_good = np.where(valid == True)
    ravalues_wrap = ravalues[index_good].copy()
    median_ra = np.nanmedian(ravalues_wrap)

    # using median to test if there is any wrapping going on
    wrap_index = np.where(np.fabs(ravalues_wrap - median_ra) > 180.0)
    nwrap = wrap_index[0].size

    # get all the ra on the same "side" of 0/360
    if(nwrap != 0 and median_ra < 180):
        ravalues_wrap[wrap_index] = ravalues_wrap[wrap_index] - 360.0

    if(nwrap != 0 and median_ra > 180):
        ravalues_wrap[wrap_index] = ravalues_wrap[wrap_index] + 360.0

    return ravalues_wrap
# ________________________________________________________________________________
