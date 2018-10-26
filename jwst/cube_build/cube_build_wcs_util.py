""" Routines related to WCS procedures of cube_build
"""
import numpy as np
from ..assign_wcs import nirspec
from gwcs import wcstools
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def find_footprint_MIRI(input, this_channel, instrument_info, coord_system):

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
       coordinate system of output cube, either alpha-beta or world

    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.
    spaxial coordinates are in units of arc seconds.
    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength

    xstart, xend = instrument_info.GetMIRISliceEndPts(this_channel)
    y, x = np.mgrid[:1024, xstart:xend]

    if coord_system == 'alpha-beta':
        detector2alpha_beta = input.meta.wcs.get_transform('detector',
                                                           'alpha_beta')
        coord1, coord2, lam = detector2alpha_beta(x, y)
    else:  # coord_system == 'world'
        coord1, coord2, lam = input.meta.wcs(x, y)
# ________________________________________________________________________________
# test for 0/360 wrapping in ra. if exists it makes it difficult to determine
# ra range of IFU cube.

    coord1_wrap = wrap_ra(coord1)
    a_min = np.nanmin(coord1_wrap)
    a_max = np.nanmax(coord1_wrap)

    b_min = np.nanmin(coord2)
    b_max = np.nanmax(coord2)

    lambda_min = np.nanmin(lam)
    lambda_max = np.nanmax(lam)

    return a_min, a_max, b_min, b_max, lambda_min, lambda_max
# ********************************************************************************


def find_footprint_NIRSPEC(input, coord_system):

    """For a NIRSPEC slice on an exposure find the foot of this data on the sky

    For each slice find:
    a. the min and max spatial coordinates (alpha,beta) or (ra,dec) depending
       on coordinate system of the output cube.
    b. min and max wavelength is also determined.

    Parameters
    ----------
    input: data model
       input model (or file)
    coord_system : str
       coordinate system of output cube, either alpha-beta or world

    Notes
    -----
    The coordinate system of alpha-beta is not yet implemented for NIRSPEC
    Returns
    -------
    min and max spaxial coordinates and wavelength for slice.

    """
    # loop over all the region (Slices) in the Channel
    # based on regions mask (indexed by slice number) find all the detector
    # x,y values for slice. Then convert the x,y values to  v2,v3,lambda
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers

    nslices = 30
    a_slice = np.zeros(nslices * 2)
    b_slice = np.zeros(nslices * 2)
    lambda_slice = np.zeros(nslices * 2)
    k = 0
    # for NIRSPEC there are 30 regions
    log.info('Looping over slices to determine cube size .. this takes a while')

    for i in range(nslices):
        slice_wcs = nirspec.nrs_wcs_set_input(input, i)
        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
        if coord_system == 'world':
            coord1, coord2, lam = slice_wcs(x, y)
        else:  # coord_system == 'alpha-beta':
            raise InvalidCoordSystem(" The Alpha-Beta Coordinate system is not valid (at this time) for NIRSPEC data")
#                detector2slicer = input.meta.wcs.get_transform('detector','slicer')
#                coord1,coord2,lam = detector2slicer(x,y)
# ________________________________________________________________________________
# For each slice  test for 0/360 wrapping in ra.
# If exists it makes it difficult to determine  ra range of IFU cube.
        coord1_wrap = wrap_ra(coord1)
        a_min = np.nanmin(coord1_wrap)
        a_max = np.nanmax(coord1_wrap)

        a_slice[k] = a_min
        a_slice[k + 1] = a_max

        b_slice[k] = np.nanmin(coord2)
        b_slice[k + 1] = np.nanmax(coord2)

        lambda_slice[k] = np.nanmin(lam)
        lambda_slice[k + 1] = np.nanmax(lam)

        k = k + 2
# ________________________________________________________________________________
# now test the ra slices for conistency. Adjust if needed.
    a_slice_wrap = wrap_ra(a_slice)
    a_min = np.nanmin(a_slice_wrap)
    a_max = np.nanmax(a_slice_wrap)

    b_min = min(b_slice)
    b_max = max(b_slice)

    lambda_min = min(lambda_slice)
    lambda_max = max(lambda_slice)

    if (a_min == 0.0 and a_max == 0.0 and b_min == 0.0 and b_max == 0.0):
        log.info('This NIRSPEC exposure has no IFU data on it - skipping file')

    return a_min, a_max, b_min, b_max, lambda_min, lambda_max
# _______________________________________________________________________________


def wrap_ra(ravalues):
    """Test for 0/360 wrapping in ra values.

    If exists it makes it difficult to determine
    ra range of IFU cube. So put them all on "one side" of 0/360 border

    Input
    -----
    ravalues : a numpy array
      ra values

    Return
    ------
    a numpy array of ra values all on "same side" of 0/360 border
    """

    valid = np.isfinite(ravalues)
    index_good = np.where(valid == True)
    ravalues_wrap = ravalues[index_good].copy()
    median_ra = np.nanmedian(ravalues_wrap)
#    print('median_ra',median_ra)

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
# Errors


class InvalidCoordSystem(Exception):
    """ Raise exeception when alpha-beta coordinate system is use for NIRSPEC
    """
    pass
