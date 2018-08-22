# Routines used for building cubes
import numpy as np
from ..assign_wcs import nirspec
from gwcs import wcstools


import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
# HELPER ROUTINES for IFUCubeData class defined in ifu_cube.py
# These methods relate to wcs type procedures.


#********************************************************************************
def find_footprint_MIRI(input, this_channel, instrument_info,coord_system):
#********************************************************************************

    """
    Short Summary
    -------------
    For each channel find:
    a. the min and max spatial coordinates (alpha,beta) or (V2-v3) depending on coordinate system.
      axis a = naxis 1, axis b = naxis2
    b. min and max wavelength is also determined. , beta and lambda for those slices


    Parameters
    ----------
    input: input model (or file)
    this_channel: channel working with


    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.
    spaxial coordinates are in units of arc seconds.
    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers

    xstart, xend = instrument_info.GetMIRISliceEndPts(this_channel)
    y, x = np.mgrid[:1024, xstart:xend]

    coord1 = np.zeros(y.shape)
    coord2 = np.zeros(y.shape)
    lam = np.zeros(y.shape)

    if coord_system == 'alpha-beta':
        detector2alpha_beta = input.meta.wcs.get_transform('detector', 'alpha_beta')
        coord1, coord2, lam = detector2alpha_beta(x, y)
    elif coord_system == 'ra-dec':
        detector2v23 = input.meta.wcs.get_transform('detector', 'v2v3')
        v23toworld = input.meta.wcs.get_transform("v2v3", "world")
        v2, v3, lam = detector2v23(x, y)
        coord1, coord2, lam = v23toworld(v2, v3, lam)
#        coord1,coord2,lam = input.meta.wcs(x,y) # for entire detector find  ra,dec,lambda
    else:
        # error the coordinate system is not defined
        raise NoCoordSystem(" The output cube coordinate system is not definded")
#________________________________________________________________________________
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

#********************************************************************************
def find_footprint_NIRSPEC(input, flag_data,coord_system):
#********************************************************************************
    """
    Short Summary
    -------------
    For each slice find:
    a. the min and max spatial coordinates (alpha,beta) or (V2-v3) depending on coordinate system.
      axis a = naxis 1, axis b = naxis2
    b. min and max wavelength is also determined. , beta and lambda for those slices


    Parameters
    ----------
    input: input model (or file)

    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.

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

        if coord_system == 'ra-dec':
            coord1, coord2, lam = slice_wcs(x, y)
        elif coord_system == 'alpha-beta':
            raise InvalidCoordSystem(" The Alpha-Beta Coordinate system is not valid (at this time) for NIRSPEC data")
#                detector2slicer = input.meta.wcs.get_transform('detector','slicer')
#                coord1,coord2,lam = detector2slicer(x,y)
        else:
            # error the coordinate system is not defined
            raise NoCoordSystem(" The output cube coordinate system is not definded")
#________________________________________________________________________________
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
#________________________________________________________________________________
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

#________________________________________________________________________________
def wrap_ra(ravalues):
    """
    Short Summary
    Test for 0/360 wrapping in ra values. If exists it makes it difficult to determine
    ra range of IFU cube. So put them all on "one side" of 0/360 border

    Input
    -----
    ravalues a numpy array of ra values

    Return
    ------
    a numpy array of ra values all on "same side" of 0/360 border
    """

    valid = np.isfinite(ravalues)
    index_good = np.where(valid == True)
#    print('number of non nan ra values',index_good[0].size,index_good[0].size/2048)
    ravalues_wrap = ravalues[index_good].copy()
    median_ra = np.nanmedian(ravalues_wrap) # find the median
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
#________________________________________________________________________________
# Errors
class NoCoordSystem(Exception):
    pass

class InvalidCoordSystem(Exception):
    pass

class RaAveError(Exception):
    pass
