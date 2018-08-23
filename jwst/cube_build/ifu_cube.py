
# Routines used for building cubes
#import sys
import time
import numpy as np
import logging
import math
from ..model_blender import blendmeta
from .. import datamodels
from ..assign_wcs import pointing
from astropy.stats import circmean
from astropy import units as u
from . import cube_build_wcs_util
from . import spaxel
from . import cube_overlap
from . import cube_cloud
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class IFUCubeData():
# IFUCubeData - holds all the important information for IFU Cube Building:
# wcs, data, reference data

    def __init__(self,
                 pipeline,
                 input_filenames,
                 input_models,
                 output_name_base,
                 output_type,
                 instrument,
                 detector,
                 list_par1,
                 list_par2,
                 instrument_info,
                 master_table,
                 **pars_cube):


        self.input_filenames = input_filenames
        self.pipeline = pipeline

        self.input_models = input_models # needed when building single mode IFU cubes
        self.output_name_base = output_name_base

        self.instrument = instrument
        self.detector = detector
        self.list_par1 = list_par1
        self.list_par2 = list_par2

        self.instrument_info = instrument_info
        self.master_table = master_table
        self.output_type = output_type
        self.scale1 = pars_cube.get('scale1')
        self.scale2 = pars_cube.get('scale2')
        self.scalew = pars_cube.get('scalew')
        self.rois = pars_cube.get('rois')
        self.roiw = pars_cube.get('roiw')

        self.interpolation = pars_cube.get('interpolation')
        self.coord_system = pars_cube.get('coord_system')
        self.offset_list = pars_cube.get('offset_list')
        self.wavemin = pars_cube.get('wavemin')
        self.wavemax = pars_cube.get('wavemax')
        self.weighting = pars_cube.get('weighting')
        self.weight_power = pars_cube.get('weight_power')
        self.xdebug = pars_cube.get('xdebug')
        self.ydebug = pars_cube.get('ydebug')
        self.zdebug = pars_cube.get('zdebug')
        self.debug_pixel = pars_cube.get('debug_pixel')
        self.spaxel_debug = pars_cube.get('spaxel_debug')

        self.num_bands = 0
        self.output_name = ''
        self.this_cube_filenames = []

        self.Cdelt1 = None
        self.Cdelt2 = None
        self.Cdelt3 = None
        self.Crpix1 = None
        self.Crpix2 = None
        self.Crpix3 = None
        self.Crval1 = None
        self.Crval2 = None
        self.Crval3 = None
        self.naxis1 = None
        self.naxis2 = None
        self.naxis3 = None

        self.a_min = 0
        self.a_max = 0
        self.b_min = 0
        self.b_max = 0
        self.lambda_min = 0
        self.lambda_max = 0
        self.xcoord = None
        self.ycoord = None
        self.zcoord = None

        self.spaxel = []        # list of spaxel classes
#********************************************************************************
    def check_ifucube(self):

        """
        Short Summary
        -------------
        Do some quick checks that the type of cube to be produced conforms to rules
        Find the ROI size.

        Parameters
        ----------
        None

        Returns
        -------
        return the ROI suze to use

        """
        num1 = len(self.list_par1)
        num_files = 0
        for i in range(num1):
            this_a = self.list_par1[i]
            this_b = self.list_par2[i]
            n = len(self.master_table.FileMap[self.instrument][this_a][this_b])
            num_files = num_files + n

# do some basic checks on the cubes
        if(self.interpolation == "area"):
            if(num_files > 1):
                raise IncorrectInput("For interpolation = area, only one file can" +
                                     " be used to created the cube")
            if(len(self.list_par1) > 1):
                raise IncorrectInput("For interpolation = area, only a single channel" +
                                     " can be used to created the cube. Use --channel=# option")
            if(self.scale2 != 0):
                raise AreaInterpolation("When using interpolation = area, the output" +
                                        " coordinate system is alpha-beta" +
                                        " The beta dimension (naxis2) has a one to one" +
                                        " mapping between slice and " +
                                        " beta coordinate.")

        if(self.coord_system == "alpha-beta"):
            if(num_files > 1):
                raise IncorrectInput("Cubes built in alpha-beta coordinate system" +
                                     " are built from a single file")
#________________________________________________________________________________
# get the ROI sizes
        roi = self.determine_roi_size()
        # if the user has not set the size of the ROI then use defaults in reference
        # parameter file

        if self.roiw == 0.0: self.roiw = roi[0]
        if self.rois == 0.0:  # user did not set so use defaults
            self.rois = roi[1]
            if self.output_type == 'single' or num_files < 4:
                self.rois = self.rois * 1.5
                log.info('Increasing spatial region of interest' +
                        'default value set for 4 dithers %f', self.rois)
        if self.interpolation == 'pointcloud':
            log.info('Region of interest spatial, wavelength  %f %f', self.rois, self.roiw)
#________________________________________________________________________________
#********************************************************************************
# update the output name
    def define_cubename(self):

#        print ('ifu_cube:define_cubename basename ',self.output_name_base)
        if self.pipeline == 2:
            newname = self.output_name_base + '_s3d.fits'
        else:
            if self.instrument == 'MIRI':
                channels = []
                for ch in self.list_par1:
                    if ch not in channels:
                        channels.append(ch)
                    number_channels = len(channels)
                    ch_name = '_ch'
                    for i in range(number_channels):
                        ch_name = ch_name + channels[i]
                        if i < number_channels - 1: ch_name = ch_name + '-'

                subchannels = list(set(self.list_par2))
                number_subchannels = len(subchannels)
                b_name = ''
                for i in range(number_subchannels):
                    b_name = b_name + subchannels[i]
                    if i > 1: b_name = b_name + '-'
                b_name = b_name.lower()
                newname = self.output_name_base + ch_name + '-' + b_name + '_s3d.fits'
                if self.coord_system == 'alpha-beta':
                    newname = self.output_name_base + ch_name + '-' + b_name + '_ab_s3d.fits'
                if self.output_type == 'single':
                    newname = self.output_name_base + ch_name + '-' + b_name + '_single_s3d.fits'
#________________________________________________________________________________
            elif self.instrument == 'NIRSPEC':
                fg_name = '_'
                for i in range(len(self.list_par1)):
                    fg_name = fg_name + self.list_par1[i] + '-' + self.list_par2[i]
                    if(i < self.num_bands - 1):
                        fg_name = fg_name + '-'
                fg_name = fg_name.lower()
                newname = self.output_name_base + fg_name + '_s3d.fits'
                if self.output_type == 'single':
                    newname = self.output_name_base + fg_name + 'single_s3d.fits'
#________________________________________________________________________________
        if self.output_type != 'single':
            log.info('Output Name %s', newname)
        return newname

#********************************************************************************
    def setup_ifucube_wcs(self):

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
        If the coordinate system is world then the min and max of
        ra(degress), dec (degrees) and lambda (microns) is returned.
        """
#________________________________________________________________________________
# Scale is 3 dimensions and is determined from values held in  instrument_info.GetScale
        scale = self.determine_scale()
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
                        ch_footprint = cube_build_wcs_util.find_footprint_NIRSPEC(
                            input_model,
                            flag_data,
                            self.coord_system)
                        amin, amax, bmin, bmax, lmin, lmax = ch_footprint
#________________________________________________________________________________
                    if self.instrument == 'MIRI':
                        ch_footprint = cube_build_wcs_util.find_footprint_MIRI(
                            input_model,
                            this_a,
                            self.instrument_info,
                            self.coord_system)
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

        if self.wavemin is not None and self.wavemin > final_lambda_min:
            final_lambda_min = self.wavemin
            log.info('Changed min wavelength of cube to %f ', final_lambda_min)

        if self.wavemax is not None and self.wavemax < final_lambda_max:
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

        if self.coord_system == 'world':
            self.set_geometry(cube_footprint)
        else:
            self.set_geometryAB(cube_footprint) # local coordinate system
        self.print_cube_geometry()

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
#_______________________________________________________________________
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
        ravalues = np.zeros(2)
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

    def set_geometryAB(self, footprint):
        """
        Short Summary
        -------------
        Based on the alpha, beta and wavelength footprint set up the size of Cube in
        alpha-beta space. This will be a single exposure cube - small FOV
        assume rectangular coord system

        """
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
                     self.naxis1, self.Crpix1, self.Crval1, self.Cdelt1,
                     self.a_min, self.a_max)
            log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                     self.naxis2, self.Crpix2, self.Crval2, self.Cdelt2,
                     self.b_min, self.b_max)
            log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                     self.naxis3, self.Crpix3, self.Crval3, self.Cdelt3,
                     self.lambda_min, self.lambda_max)

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

#********************************************************************************
    def build_ifucube(self):

        """
        Short Summary
        -------------
        Loop over every band contained in the IFU cube and read in the data associated with the band
        Map the detector data to the cube output coordinate system

        Parameter
        ----------
        spaxel - a list of spaxel members holding the detector flux information

        Returns
        -------
        each spaxel element with the mapped detector values associated with it

        """

        self.output_name = self.define_cubename()
        self.spaxel = self.create_spaxel()

        # now need to loop over every file that covers this channel/subchannel (MIRI)
        # or Grating/filter(NIRSPEC)
        #and map the detector pixels to the cube spaxel.

        number_bands = len(self.list_par1)
        t0 = time.time()
        for i in range(number_bands):
            this_par1 = self.list_par1[i]
            this_par2 = self.list_par2[i]

            log.debug("Working on Band defined by:%s %s ", this_par1, this_par2)
            self.map_detector_to_spaxel(this_par1, this_par2, self.spaxel)

        t1 = time.time()
        log.info("Time Map All slices on Detector to Cube = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
# Mapped all data to cube or Point Cloud
# now determine Cube Spaxel flux

        t0 = time.time()
        self.find_spaxel_flux(self.spaxel)

        t1 = time.time()
        log.info("Time to find Cube Flux= %.1f.s" % (t1 - t0,))

        ifucube_model = self.setup_ifucube(0)
#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube
        self.update_ifucube(ifucube_model, self.spaxel)
        return ifucube_model

#********************************************************************************

    def build_ifucube_single(self):

        """
        Short Summary
        -------------
        Loop over every band contained in the IFU cube and read in the data associated with the band
        Map the detector data to the cube output coordinate system

        Parameter
        ----------
        spaxel - a list of spaxel members holding the detector flux information

        Returns
        -------
        each spaxel element with the mapped detector values associated with it

        """

        # loop over input models

        single_ifucube_container = datamodels.ModelContainer()
        n = len(self.input_models)
        log.info("Number of Single IFU cubes creating  = %i" % n)
        this_par1 = self.list_par1[0] # only one channel is used in this approach
        this_par2 = None # not important for this type of mapping

        self.weighting == 'msm'
        c1_offset = 0
        c2_offset = 0
        for j in range(n):
            log.info("Working on next Single IFU Cube  = %i" % (j + 1))
            t0 = time.time()
# for each new data model create a new spaxel
            spaxel = []
            spaxel = self.create_spaxel()
            with datamodels.IFUImageModel(self.input_models[j]) as input_model:
#********************************************************************************
# pulled necessary routines from   CubeData.map_detector_to_spaxel
                if self.instrument == 'MIRI':
#________________________________________________________________________________
                    xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_par1)
                    y, x = np.mgrid[:1024, xstart:xend]
                    y = np.reshape(y, y.size)
                    x = np.reshape(x, x.size)

                    cube_cloud.match_det2cube(self, input_model,
                                              x, y, j,
                                              this_par1, this_par2,
                                              spaxel,
                                              c1_offset, c2_offset)

                elif self.instrument == 'NIRSPEC':
                    # each file, detector has 30 slices - wcs information access seperately for each slice
                    nslices = 30
                    for ii in range(nslices):
                        t0a = time.time()
                        #slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)
                        #x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
                        x = None
                        y = None

                        cube_cloud.match_det2cube(self, input_model,
                                                  x, y, ii,
                                                  this_par1, this_par2,
                                                  spaxel,
                                                  c1_offset, c2_offset)

                        t1a = time.time()
                        log.debug("Time Match one NIRSPEC slice to ifucube = %.1f.s" % (t1a - t0a,))
#_______________________________________________________________________
# shove Flux and iflux in the  final ifucube

            self.find_spaxel_flux(spaxel)
# now determine Cube Spaxel flux

            ifucube_model = self.setup_ifucube(j)
            self.update_ifucube(ifucube_model, spaxel)
            t1 = time.time()
            log.info("Time Create Single ifucube  = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
            single_ifucube_container.append(ifucube_model)
            del spaxel[:]
        return single_ifucube_container
#********************************************************************************
    def create_spaxel(self):
        """
        Short Summary
        -------------
        # now you have the size of cube - create an instance for each spaxel
        # create an empty spaxel list - this will become a list of Spaxel classses

        Parameter
        ----------

        Returns
        -------
        list of classes contained in spaxel
        """
#________________________________________________________________________________
        total_num = self.naxis1 * self.naxis2 * self.naxis3

        if(self.interpolation == 'pointcloud'):
            for t in range(total_num):
                self.spaxel.append(spaxel.Spaxel())
        else:
            for t in range(total_num):
                self.spaxel.append(spaxel.SpaxelAB())

        return self.spaxel

#********************************************************************************
    def determine_roi_size(self):
        """
        Short Summary
        -------------
        Determine the spatial and wavelength roi size to use for selecting point cloud
        elements around the spaxel centeres.
        If the IFU cube covers more than 1 band - then use the rules to
        define the Spatial and Wavelength roi size to use for the cube
        Current Rule: using the minimum

        Parameters
        ----------
        instrument_info holds the defaults roi sizes  for each channel/subchannel (MIRI)
        or grating (NIRSPEC)

        Returns
        -------
        roi size for spatial and wavelength

        """

        roi = [0, 0]
        if self.instrument == 'MIRI':
            number_bands = len(self.list_par1)
            min_s = 1000.00
            min_w = 1000.00

            for i in range(number_bands):
                this_channel = self.list_par1[i]
                this_sub = self.list_par2[i]
                wroi = self.instrument_info.GetWaveRoi(this_channel, this_sub)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_channel, this_sub)
                if sroi < min_s:
                    min_s = sroi
            roi = [min_w, min_s]

        elif self.instrument == 'NIRSPEC':
            number_gratings = len(self.list_par1)
            min_s = 1000.00
            min_w = 1000.00

            for i in range(number_gratings):
                this_gwa = self.list_par1[i]
                this_filter = self.list_par2[i]

                wroi = self.instrument_info.GetWaveRoi(this_gwa, this_filter)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_gwa, this_filter)
                if sroi < min_s:
                    min_s = sroi
            roi = [min_w, min_s]
        return roi

#********************************************************************************
    def map_detector_to_spaxel(self, this_par1, this_par2, spaxel):
        from ..mrs_imatch.mrs_imatch_step import apply_background_2d
#********************************************************************************
        """
        Short Summary
        -------------
        Loop over files that cover the cube and map the detector pixel to Cube spaxels
        If dither offsets have been supplied then apply those values to the data

        Parameter
        ----------
        spaxel: List of Spaxels

        Returns
        -------
        if(interpolation = area - only valid for alpha-beta
        or
        if(interpolation = pointcloud
        """

        instrument = self.instrument
        nfiles = len(self.master_table.FileMap[instrument][this_par1][this_par2])

    # loop over the files that cover the spectral range the cube is for

        for k in range(nfiles):
            ifile = self.master_table.FileMap[instrument][this_par1][this_par2][k]
            ioffset = len(self.master_table.FileOffset[this_par1][this_par2]['C1'])

            self.this_cube_filenames.append(ifile)

            c1_offset = 0.0
            c2_offset = 0.0
        # c1_offset and c2_offset are the dither offset sets (in arc seconds)
        # by default these are zer0. The user has to supply these
            if ioffset == nfiles:
                c1_offset = self.master_table.FileOffset[this_par1][this_par2]['C1'][k]
                c2_offset = self.master_table.FileOffset[this_par1][this_par2]['C2'][k]
# Open the input data model
            with datamodels.IFUImageModel(ifile) as input_model:
                # check if background sky matching as been done
                # mrs_imatch step. THis is only for MRS data at this time
                # but go head and check it before splitting by instrument
                # the polynomial should be empty for NIRSPEC
                #do_background_subtraction = False
                num_ch_bgk = len(input_model.meta.background.polynomial_info)

                if(num_ch_bgk > 0):

                    #do_background_subtraction = True
                    for ich_num in range(num_ch_bgk):
                        poly = input_model.meta.background.polynomial_info[ich_num]
                        poly_ch = poly.channel
                        if(poly_ch == this_par1):
                            apply_background_2d(input_model, poly_ch, subtract=True)

#********************************************************************************
                if self.instrument == 'MIRI':
#________________________________________________________________________________
# Standard method
                    if self.interpolation == 'pointcloud':
                        xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_par1)
                        y, x = np.mgrid[:1024, xstart:xend]
                        y = np.reshape(y, y.size)
                        x = np.reshape(x, x.size)
                        t0 = time.time()
                        cube_cloud.match_det2cube(self, input_model,
                                            x, y, k,
                                            this_par1, this_par2,
                                            spaxel,
                                            c1_offset, c2_offset)


                        t1 = time.time()
                        log.debug("Time Match one Channel from 1 file to ifucube = %.1f.s"
                                  % (t1 - t0,))
#________________________________________________________________________________
#2D area method - only works for single files and coord_system = 'alpha-beta'
                    if self.interpolation == 'area':
                        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                              'alpha_beta')


                        start_region = self.instrument_info.GetStartSlice(this_par1)
                        end_region = self.instrument_info.GetEndSlice(this_par1)
                        regions = list(range(start_region, end_region + 1))

                        for i in regions:
                            log.info('Working on Slice # %d', i)

                            y, x = (det2ab_transform.label_mapper.mapper == i).nonzero()

                    # spaxel object holds all needed information in a set of lists
                    #    flux (of overlapping detector pixel)
                    #    error (of overlapping detector pixel)
                    #    overlap ratio
                    #    beta distance

# getting pixel corner - ytop = y + 1 (routine fails for y = 1024)
                            index = np.where(y < 1023)
                            y = y[index]
                            x = x[index]
                            t0 = time.time()

                            cube_overlap.match_det2cube(self, x, y, i,
                                                        start_region,
                                                        input_model,
                                                        det2ab_transform,
                                                        spaxel)
                            t1 = time.time()
                            log.debug("Time Map one Slice  to Cube = %.1f.s" % (t1 - t0,))

#********************************************************************************
                elif instrument == 'NIRSPEC':
                    # each file, detector has 30 slices - wcs information access seperately for each slice
                    nslices = 30
                    log.info("Mapping each NIRSPEC slice to sky, this takes a while for NIRSPEC data")
                    for i in range(nslices):
#                        slice_wcs = nirspec.nrs_wcs_set_input(input_model, i)
#                        x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box,
#                                                              step=(1,1), center=True)

                        t0 = time.time()
                        x = 0
                        y = 0
                        cube_cloud.match_det2cube(self, input_model,
                                                  x, y, i,
                                                  this_par1, this_par2,
                                                  spaxel,
                                                  c1_offset, c2_offset)


                        t1 = time.time()
                        log.debug("Time Match one NIRSPEC slice to ifucube = %.1f.s" % (t1 - t0,))
#********************************************************************************
    def find_spaxel_flux(self, spaxel):
#********************************************************************************
        """
        Short Summary
        -------------
        Depending on the interpolation method, find the flux for each spaxel value

        Parameter
        ----------
        spaxel: List of Spaxels
        PixelCloud - pixel point cloud, only filled in if doing 3-D interpolation

        Returns
        -------
        if(interpolation = area) flux determined for each spaxel
        or
        if(interpolation = pointcloud) flux determined for each spaxel based on interpolation of PixelCloud
        """


        if self.interpolation == 'area':
            nspaxel = len(spaxel)

            for i in range(nspaxel):
                if(spaxel[i].iflux > 0):
                    spaxel[i].flux = spaxel[i].flux / spaxel[i].flux_weight

        elif self.interpolation == 'pointcloud':
            icube = 0
            t0 = time.time()
            for iz, z in enumerate(self.zcoord):
                for iy, y in enumerate(self.ycoord):
                    for ix, x in enumerate(self.xcoord):

                        if(spaxel[icube].iflux > 0):
                            spaxel[icube].flux = spaxel[icube].flux / spaxel[icube].flux_weight

#                            if(self.debug_pixel == 1 and self.xdebug == ix and
#                               self.ydebug == iy and self.zdebug == iz ):
#                                log.debug('For spaxel %d %d %d final flux %f '
#                                          %(self.xdebug + 1,self.ydebug + 1,
#                                            self.zdebug + 1,spaxel[icube].flux))
#                                self.spaxel_debug.write('For spaxel %d %d %d, final flux %f '
#                                                        %(self.xdebug + 1,self.ydebug + 1,
#                                                          self.zdebug + 1,spaxel[icube].flux) +' \n')
                        icube = icube + 1
            t1 = time.time()
            log.info("Time to interpolate at spaxel values = %.1f.s" % (t1 - t0,))

#********************************************************************************
    def setup_ifucube(self, j):

        """
        Short Summary
        -------------
        Set up the final  the IFU cube to fits file

        Parameters
        ----------
        Cube: holds meta data of cube
        spaxel: list of spaxels in cube


        Returns
        -------
        return IFUCube model

        """
        naxis1 = self.naxis1
        naxis2 = self.naxis2
        naxis3 = self.naxis3

        data = np.zeros((naxis3, naxis2, naxis1))
        idata = np.zeros((naxis3, naxis2, naxis1))

        dq_cube = np.zeros((naxis3, naxis2, naxis1))
        err_cube = np.zeros((naxis3, naxis2, naxis1))

        ifucube_model = datamodels.IFUCubeModel(data=data, dq=dq_cube, err=err_cube, weightmap=idata)
        ifucube_model.update(self.input_models[j])
        ifucube_model.meta.filename = self.output_name

        # Call model_blender if there are multiple inputs
        if len(self.input_models) > 1:
            saved_model_type = ifucube_model.meta.model_type
            self.blend_output_metadata(ifucube_model)
            ifucube_model.meta.model_type = saved_model_type  # Reset to original

#______________________________________________________________________
        if self.output_type == 'single':
            with datamodels.open(self.input_models[j]) as input:
                # define the cubename for each single
                filename = self.input_filenames[j]
                indx = filename.rfind('.fits')
                self.output_name_base = filename[:indx]
                self.output_file = None
                newname = self.define_cubename()
                ifucube_model.meta.filename = newname

#______________________________________________________________________
# fill in Channel, Band for MIRI
        if self.instrument == 'MIRI':
            # fill in Channel output meta data
            num_ch = len(self.list_par1)
            ifucube_model.meta.instrument.channel = self.list_par1[0]
            num_ch = len(self.list_par1)
            for m in range(1, num_ch):
                ifucube_model.meta.instrument.channel = (
                ifucube_model.meta.instrument.channel + str(self.list_par1[m]))

#______________________________________________________________________
        ifucube_model.meta.wcsinfo.crval1 = self.Crval1
        ifucube_model.meta.wcsinfo.crval2 = self.Crval2
        ifucube_model.meta.wcsinfo.crval3 = self.Crval3
        ifucube_model.meta.wcsinfo.crpix1 = self.Crpix1
        ifucube_model.meta.wcsinfo.crpix2 = self.Crpix2
        ifucube_model.meta.wcsinfo.crpix3 = self.Crpix3
        ifucube_model.meta.wcsinfo.cdelt1 = self.Cdelt1 / 3600.0
        ifucube_model.meta.wcsinfo.cdelt2 = self.Cdelt2 / 3600.0
        ifucube_model.meta.wcsinfo.cdelt3 = self.Cdelt3

        ifucube_model.meta.wcsinfo.ctype1 = 'RA---TAN'
        ifucube_model.meta.wcsinfo.ctype2 = 'DEC--TAN'
        ifucube_model.meta.wcsinfo.cunit1 = 'deg'
        ifucube_model.meta.wcsinfo.cunit2 = 'deg'

        ifucube_model.meta.wcsinfo.ctype3 = 'WAVE'
        ifucube_model.meta.wcsinfo.cunit3 = 'um'
        ifucube_model.meta.wcsinfo.wcsaxes = 3
        ifucube_model.meta.wcsinfo.pc1_1 = -1
        ifucube_model.meta.wcsinfo.pc1_2 = 0
        ifucube_model.meta.wcsinfo.pc1_3 = 0

        ifucube_model.meta.wcsinfo.pc2_1 = 0
        ifucube_model.meta.wcsinfo.pc2_2 = 1
        ifucube_model.meta.wcsinfo.pc2_3 = 0

        ifucube_model.meta.wcsinfo.pc3_1 = 0
        ifucube_model.meta.wcsinfo.pc3_2 = 0
        ifucube_model.meta.wcsinfo.pc3_3 = 1

        ifucube_model.meta.ifu.flux_extension = 'SCI'
        ifucube_model.meta.ifu.error_extension = 'ERR'
        ifucube_model.meta.ifu.error_type = 'ERR'
        ifucube_model.meta.ifu.dq_extension = 'DQ'
        ifucube_model.meta.ifu.roi_spatial = self.rois
        ifucube_model.meta.ifu.roi_wave = self.roiw
        ifucube_model.meta.ifu.weighting = self.weighting
        ifucube_model.meta.ifu.weight_power = self.weight_power

        with datamodels.open(self.input_models[j]) as input:
            ifucube_model.meta.bunit_data = input.meta.bunit_data
            ifucube_model.meta.bunit_err = input.meta.bunit_err

        if self.coord_system == 'alpha-beta':
            ifucube_model.meta.wcsinfo.cunit1 = 'arcsec'
            ifucube_model.meta.wcsinfo.cunit2 = 'arcsec'

# we only need to check list_par1[0] and list_par2[0] because these types
# of cubes are made from 1 exposures (setup_cube checks this at the start
# of cube_build).
            if self.list_par1[0] == '1' and self.list_par2[0] == 'short':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL1A'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE1A'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'short':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL2A'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE2A'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'short':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL3A'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE3A'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'short':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL4A'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE4A'

            if self.list_par1[0] == '1' and self.list_par2[0] == 'medium':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL1B'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE1B'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'medium':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL2B'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE2B'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'medium':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL3B'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE3B'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'medium':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL4B'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE4B'

            if self.list_par1[0] == '1' and self.list_par2[0] == 'long':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL1C'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE1C'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'long':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL2C'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE2C'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'long':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL3C'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE3C'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'long':
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSAL4C'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBE4C'


        wcsobj = pointing.create_fitswcs(ifucube_model)
        ifucube_model.meta.wcs = wcsobj
        ifucube_model.meta.wcs.bounding_box = ((0, naxis1 - 1),
                                         (0, naxis2 - 1),
                                         (0, naxis3 - 1))
        return ifucube_model

#********************************************************************************
    def update_ifucube(self, ifucube_model, spaxel):
#********************************************************************************
        """
        Short Summary
        -------------
        Write the IFU cube to fits file

        Parameters
        ----------
        Cube: holds meta data of cube
        spaxel: list of spaxels in cube


        Returns
        -------
        fills in IFUdata arrays with spaxel

        """
    #pull out data into array
        temp_flux = np.reshape(np.array([s.flux for s in spaxel]),
                          [self.naxis3, self.naxis2, self.naxis1])
        temp_wmap = np.reshape(np.array([s.iflux for s in spaxel]),
                          [self.naxis3, self.naxis2, self.naxis1])


        ifucube_model.data = temp_flux
        ifucube_model.weightmap = temp_wmap
        ifucube_model.meta.cal_step.cube_build = 'COMPLETE'
#    icube = 0
#    for z in range(Cube.naxis3):
#        for y in range(Cube.naxis2):
#            for x in range(Cube.naxis1):
#                IFUCube.data[z, y, x] = spaxel[icube].flux
#                IFUCube.weightmap[z, y, x] = len(spaxel[icube].ipointcloud)
#                icube = icube + 1

#********************************************************************************
    def blend_output_metadata(self, IFUCube):

        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = IFUCube.meta.filename
#        log.info('Blending metadata for {}'.format(output_file))
        blendmeta.blendmodels(IFUCube, inputs=self.input_models,
                              output=output_file)
class IncorrectInput(Exception):
    pass

class AreaInterpolation(Exception):
    pass
