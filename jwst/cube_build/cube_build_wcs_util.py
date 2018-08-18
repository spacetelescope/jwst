# Routines used for building cubes
import numpy as np
import math
from .. import datamodels
from ..assign_wcs import nirspec
from . import coord
from gwcs import wcstools
from astropy.stats import circmean
from astropy import units as u

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
# HELPER ROUTINES for IFUCubeData class defined in ifu_cube.py
# These methods relate to wcs type procedures.

#********************************************************************************
def setup_wcs(self):
#********************************************************************************
    """
    Short Summary
    -------------
    Function to determine the min and max coordinates of the spectral
    cube, given channel & subchannel


    Parameter
    ----------
    self.master_table:  A table that contains the channel/subchannel or
    filter/grating for each input file
    self.instrument_info: Default information on the MIRI and NIRSPEC instruments.

    Returns
    -------
    Cube Dimension Information:
    Footprint of cube: min and max of coordinates of cube.
    If an offset list is provided then these values are applied.
    If the coordinate system is alpha-beta (MIRI) then min and max
    coordinates of alpha (arc sec), beta (arc sec) and lambda (microns)
    If the coordinate system is ra-dec then the min and max of
    ra(degress), dec (degrees) and lambda (microns) is returned.
    """

#________________________________________________________________________________
# Scale is 3 dimensions and is determined from values held in  instrument_info.GetScale
    scale = determine_scale(self)
    self.Cdelt1 = scale[0]
    self.Cdelt2 = scale[1]
    self.Cdelt3 = scale[2]

    parameter1 = self.list_par1
    parameter2 = self.list_par2
    a_min = []
    a_max = []
    b_min = []
    b_max = []
    lambda_min = []
    lambda_max = []

    self.num_bands = len(self.list_par1)
    log.info('Number of bands in cube  %i', self.num_bands)

    for i in range(self.num_bands):
        this_a = parameter1[i]
        this_b = parameter2[i]
        log.debug('Working on data  from %s,%s', this_a, this_b)
        n = len(self.master_table.FileMap[self.instrument][this_a][this_b])
        log.debug('number of files %d ', n)
    # each file find the min and max a and lambda (OFFSETS NEED TO BE APPLIED TO THESE VALUES)
        for k in range(n):
            amin = 0.0
            amax = 0.0
            bmin = 0.0
            bmax = 0.0
            lmin = 0.0
            lmax = 0.0
            c1_offset = 0.0
            c2_offset = 0.0
            ifile = self.master_table.FileMap[self.instrument][this_a][this_b][k]
            ioffset = len(self.master_table.FileOffset[this_a][this_b]['C1'])
            if ioffset == n:
                c1_offset = self.master_table.FileOffset[this_a][this_b]['C1'][k]
                c2_offset = self.master_table.FileOffset[this_a][this_b]['C2'][k]
#________________________________________________________________________________
# Open the input data model
# Find the footprint of the image

            with datamodels.IFUImageModel(ifile) as input_model:
                if self.instrument == 'NIRSPEC':
                    flag_data = 0
                    ch_footprint = find_footprint_NIRSPEC(self,
                                                          input_model,
                                                          flag_data)
                    amin, amax, bmin, bmax, lmin, lmax = ch_footprint
#________________________________________________________________________________
                if self.instrument == 'MIRI':
                    ch_footprint = find_footprint_MIRI(self,
                                                       input_model,
                                                       this_a,
                                                       self.instrument_info)
                    amin, amax, bmin, bmax, lmin, lmax = ch_footprint

# If a dither offset list exists then apply the dither offsets (offsets in arc seconds)

                amin = amin - c1_offset / 3600.0
                amax = amax - c1_offset / 3600.0

                bmin = bmin - c2_offset / 3600.0
                bmax = bmax - c2_offset / 3600.0

                a_min.append(amin)
                a_max.append(amax)
                b_min.append(bmin)
                b_max.append(bmax)
                lambda_min.append(lmin)
                lambda_max.append(lmax)
#________________________________________________________________________________
    # done looping over files determine final size of cube

    final_a_min = min(a_min)
    final_a_max = max(a_max)
    final_b_min = min(b_min)
    final_b_max = max(b_max)
    final_lambda_min = min(lambda_min)
    final_lambda_max = max(lambda_max)

    if self.wavemin != None and self.wavemin > final_lambda_min:
        final_lambda_min = self.wavemin
        log.info('Changed min wavelength of cube to %f ', final_lambda_min)

    if self.wavemax != None and self.wavemax < final_lambda_max:
        final_lambda_max = self.wavemax
        log.info('Changed max wavelength of cube to %f ', final_lambda_max)
#________________________________________________________________________________
    if self.instrument == 'MIRI' and self.coord_system == 'alpha-beta':
        #  we have a 1 to 1 mapping in beta dimension.
        nslice = self.instrument_info.GetNSlice(parameter1[0])
        log.info('Beta Scale %f ', self.Cdelt2)
        self.Cdelt2 = (final_b_max - final_b_min) / nslice
        final_b_max = final_b_min + (nslice) * self.Cdelt2
        log.info('Changed the Beta Scale dimension so we have 1 -1 mapping between beta and slice #')
        log.info('New Beta Scale %f ', self.Cdelt2)
#________________________________________________________________________________
# Test that we have data (NIRSPEC NRS2 only has IFU data for 3 configurations)
    test_a = final_a_max - final_a_min
    test_b = final_b_max - final_b_min
    test_w = final_lambda_max - final_lambda_min
    tolerance1 = 0.00001
    tolerance2 = 0.1
    if(test_a < tolerance1 or test_b < tolerance1 or test_w < tolerance2):
        log.info('No Valid IFU slice data found %f %f %f ', test_a, test_b, test_w)
#________________________________________________________________________________
    cube_footprint = (final_a_min, final_a_max, final_b_min, final_b_max,
                      final_lambda_min, final_lambda_max)
#________________________________________________________________________________
    # Based on Scaling and Min and Max values determine naxis1, naxis2, naxis3
    # set cube CRVALs, CRPIXs and xyz coords (center  x,y,z vector spaxel centers)

    if self.coord_system == 'ra-dec':
        set_geometry(self, cube_footprint)
    else:
        set_geometryAB(self, cube_footprint) # local coordinate system
    print_cube_geometry(self)

#********************************************************************************
def determine_scale(self):
#********************************************************************************
    """
    Short Summary
    -------------
    Determine the scale (sampling) in the 3 dimensions for the cube.
    If the IFU cube covers more than 1 band - then use the rules to
    define the spatial and spectral sample size to use for the cube
    Current Rule: using the minimum

    Parameters
    ----------
    self.instrument_info holds the defaults scales for each channel/subchannel (MIRI)
    or Grating (NIRSPEC)

    Returns
    -------
    scale, holding the scale for the 3 dimensions of the cube/

    """
    scale = [0, 0, 0]
    if self.instrument == 'MIRI':
        number_bands = len(self.list_par1)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_bands):
            this_channel = self.list_par1[i]
            this_sub = self.list_par2[i]
            a_scale, b_scale, w_scale = self.instrument_info.GetScale(this_channel,
                                                                      this_sub)
            if a_scale < min_a:
                min_a = a_scale
            if b_scale < min_b:
                min_b = b_scale
            if w_scale < min_w:
                min_w = w_scale
        scale = [min_a, min_b, min_w]

    elif self.instrument == 'NIRSPEC':
        number_gratings = len(self.list_par1)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_gratings):
            this_gwa = self.list_par1[i]
            this_filter = self.list_par2[i]
            a_scale, b_scale, w_scale = self.instrument_info.GetScale(this_gwa,
                                                                      this_filter)
            if a_scale < min_a:
                min_a = a_scale
            if b_scale < min_b:
                min_b = b_scale
            if w_scale < min_w:
                min_w = w_scale
        scale = [min_a, min_b, min_w]

#________________________________________________________________________________
# check and see if the user has set the scale or set by cfg.

    a_scale = scale[0]
    if self.scale1 != 0.0:
        a_scale = self.scale1

    b_scale = scale[1]
    if self.scale2 != 0.0:
        b_scale = self.scale2

    w_scale = scale[2]
        # temp fix for large cubes - need to change to variable wavelength scale
    if self.scalew == 0 and self.num_bands > 6:
        w_scale = w_scale * 2
    if self.scalew == 0 and self.num_bands > 9:
        w_scale = w_scale * 2
    if self.scalew != 0.0:
        w_scale = self.scalew

    scale = [a_scale, b_scale, w_scale]
    return scale
#_______________________________________________________________________


#********************************************************************************
def find_footprint_MIRI(self, input, this_channel, instrument_info):
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

    if self.coord_system == 'alpha-beta':
        detector2alpha_beta = input.meta.wcs.get_transform('detector', 'alpha_beta')
        coord1, coord2, lam = detector2alpha_beta(x, y)
    elif self.coord_system == 'ra-dec':
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
def find_footprint_NIRSPEC(self, input, flag_data):
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

        if self.coord_system == 'ra-dec':
            coord1, coord2, lam = slice_wcs(x, y)

        elif self.coord_system == 'alpha-beta':
            raise InvalidCoordSystem(" The Alpha-Beta Coordinate system is not valid (at this time) for NIRSPEC data")
#                detector2slicer = input.meta.wcs.get_transform('detector','slicer')
#                coord1,coord2,lam = detector2slicer(x,y)

        else:
            # error the coordinate system is not defined
            raise NoCoordSystem(" The output cube coordinate system is not definded")
#________________________________________________________________________________
# For each slice  test for 0/360 wrapping in ra. 
# If exists it makes it difficult to determine  ra range of IFU cube. 
##            print(' # ra values',ra.size,ra.size/2048)
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
#_______________________________________________________________________
# Footprint values are RA,DEC values on the sky
# Values are given in degrees

def set_geometry(self, footprint):

    """
    Short Summary
    -------------
    Based on the ra,dec and wavelength footprint set up the size of Cube in
    the tangent plane projected coordinate system.

    """

    ra_min, ra_max, dec_min, dec_max, lambda_min, lambda_max = footprint # in degrees
    dec_ave = (dec_min + dec_max) / 2.0

        # we can not average ra values because of the convergence of hour angles.
    ravalues = np.zeros(2) # we might want to increase the number of ravalues later
                                # is just taking min and max is not sufficient
    ravalues[0] = ra_min
    ravalues[1] = ra_max

        # astropy circmean assumes angles are in radians, we have angles in degrees
    ra_ave = circmean(ravalues * u.deg).value
    log.info('Ra average %f12.8', ra_ave)

    self.Crval1 = ra_ave
    self.Crval2 = dec_ave
    xi_center, eta_center = coord.radec2std(self.Crval1, self.Crval2,
                                           ra_ave, dec_ave)
    xi_min, eta_min = coord.radec2std(self.Crval1, self.Crval2, ra_min, dec_min)
    xi_max, eta_max = coord.radec2std(self.Crval1, self.Crval2, ra_max, dec_max)

#________________________________________________________________________________
# find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
# to find location of center abs of min values is how many pixels

    n1a = int(math.ceil(math.fabs(xi_min) / self.Cdelt1))
    n2a = int(math.ceil(math.fabs(eta_min) / self.Cdelt2))

    n1b = int(math.ceil(math.fabs(xi_max) / self.Cdelt1))
    n2b = int(math.ceil(math.fabs(eta_max) / self.Cdelt2))

    xi_min = 0.0 - (n1a * self.Cdelt1) - (self.Cdelt1 / 2.0)
    xi_max = (n1b * self.Cdelt1) + (self.Cdelt1 / 2.0)

    eta_min = 0.0 - (n2a * self.Cdelt2) - (self.Cdelt2 / 2.0)
    eta_max = (n2b * self.Cdelt2) + (self.Cdelt2 / 2.0)

    self.Crpix1 = float(n1a) + 1.0
    self.Crpix2 = float(n2a) + 1.0

    self.naxis1 = n1a + n1b
    self.naxis2 = n2a + n2b

    self.a_min = xi_min
    self.a_max = xi_max
    self.b_min = eta_min
    self.b_max = eta_max

# center of spaxels
    self.xcoord = np.zeros(self.naxis1)
    xstart = xi_min + self.Cdelt1 / 2.0
    for i in range(self.naxis1):
        self.xcoord[i] = xstart
        xstart = xstart + self.Cdelt1

    self.ycoord = np.zeros(self.naxis2)
    ystart = eta_min + self.Cdelt2 / 2.0

    for i in range(self.naxis2):
        self.ycoord[i] = ystart
        ystart = ystart + self.Cdelt2

    ygrid = np.zeros(self.naxis2 * self.naxis1)
    xgrid = np.zeros(self.naxis2 * self.naxis1)

    k = 0
    ystart = self.ycoord[0]
    for i in range(self.naxis2):
        xstart = self.xcoord[0]
        for j in range(self.naxis1):
            xgrid[k] = xstart
            ygrid[k] = ystart
            xstart = xstart + self.Cdelt1
            k = k + 1
        ystart = ystart + self.Cdelt2


    self.Xcenters = xgrid
    self.Ycenters = ygrid
#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

    self.lambda_min = lambda_min
    self.lambda_max = lambda_max
    range_lambda = self.lambda_max - self.lambda_min
    self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
    lambda_center = (self.lambda_max + self.lambda_min) / 2.0
    self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
    self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

    self.zcoord = np.zeros(self.naxis3)
    self.Crval3 = self.lambda_min
    self.Crpix3 = 1.0
    zstart = self.lambda_min + self.Cdelt3 / 2.0

    for i in range(self.naxis3):
        self.zcoord[i] = zstart
        zstart = zstart + self.Cdelt3
#_______________________________________________________________________
# cube in alpha-beta space (single exposure cube - small FOV assume rectangular coord system
def set_geometryAB(self, footprint):
    self.a_min, self.a_max, self.b_min, self.b_max, self.lambda_min, self.lambda_max = footprint

        #set up the a (x) coordinates of the cube
    range_a = self.a_max - self.a_min
    self.naxis1 = int(math.ceil(range_a / self.Cdelt1))

        # adjust min and max based on integer value of naxis1
    a_center = (self.a_max + self.a_min) / 2.0
    self.a_min = a_center - (self.naxis1 / 2.0) * self.Cdelt1
    self.a_max = a_center + (self.naxis1 / 2.0) * self.Cdelt1

    self.xcoord = np.zeros(self.naxis1)
    self.Crval1 = self.a_min
    self.Crpix1 = 0.5
    xstart = self.a_min + self.Cdelt1 / 2.0
    for i in range(self.naxis1):
        self.xcoord[i] = xstart
        xstart = xstart + self.Cdelt1

#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube

    range_lambda = self.lambda_max - self.lambda_min
    self.naxis3 = int(math.ceil(range_lambda / self.Cdelt3))

         # adjust max based on integer value of naxis3
    lambda_center = (self.lambda_max + self.lambda_min) / 2.0

    self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.Cdelt3
    self.lambda_max = lambda_center + (self.naxis3 / 2.0) * self.Cdelt3

    self.lambda_max = self.lambda_min + (self.naxis3) * self.Cdelt3

    self.zcoord = np.zeros(self.naxis3)
    self.Crval3 = self.lambda_min
    self.Crpix3 = 1.0
    zstart = self.lambda_min + self.Cdelt3 / 2.0

    for i in range(self.naxis3):
        self.zcoord[i] = zstart
        zstart = zstart + self.Cdelt3
#_______________________________________________________________________
        # set up the naxis2 parameters
    range_b = self.b_max - self.b_min

    self.naxis2 = int(math.ceil(range_b / self.Cdelt2))
    b_center = (self.b_max + self.b_min) / 2.0

        # adjust min and max based on integer value of naxis2
    self.b_max = b_center + (self.naxis2 / 2.0) * self.Cdelt2
    self.b_min = b_center - (self.naxis2 / 2.0) * self.Cdelt2

    self.ycoord = np.zeros(self.naxis2)
    self.Crval2 = self.b_min
    self.Crpix2 = 0.5
    ystart = self.b_min + self.Cdelt2 / 2.0
    for i in range(self.naxis2):
        self.ycoord[i] = ystart
        ystart = ystart + self.Cdelt2

#_______________________________________________________________________

def print_cube_geometry(self):

    """
    Print out the general properties of the size of the IFU Cube
    """

    log.info('Cube Geometry:')
    if self.coord_system == 'alpha-beta':
        log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (alpha,beta arc sec)')
    else:
        log.info('axis# Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (xi,eta arc sec)')
    log.info('Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
             self.naxis1, self.Crpix1, self.Crval1, self.Cdelt1, self.a_min, self.a_max)
    log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
             self.naxis2, self.Crpix2, self.Crval2, self.Cdelt2, self.b_min, self.b_max)
    log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
             self.naxis3, self.Crpix3, self.Crval3, self.Cdelt3, self.lambda_min, self.lambda_max)

    if self.instrument == 'MIRI':
        # length of channel and subchannel are the same
        number_bands = len(self.list_par1)

        for i in range(number_bands):
            this_channel = self.list_par1[i]
            this_subchannel = self.list_par2[i]
            log.info('Cube covers channel, subchannel: %s %s ', this_channel, this_subchannel)
    elif self.instrument == 'NIRSPEC':
            # number of filters and gratings are the same
        number_bands = len(self.list_par1)

        for i in range(number_bands):
            this_fwa = self.list_par2[i]
            this_gwa = self.list_par1[i]
            log.info('Cube covers grating, filter: %s %s ', this_gwa, this_fwa)
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

    index_good = np.where( valid == True)
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
