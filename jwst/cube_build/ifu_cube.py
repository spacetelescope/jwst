
# Routines used for building cubes
import time
import numpy as np
import logging
import math
from ..model_blender import blendmeta
from .. import datamodels
from ..assign_wcs import pointing
from jwst.transforms.models import _toindex
from astropy.stats import circmean
from astropy import units as u
from gwcs import wcstools
from ..assign_wcs import nirspec
from ..datamodels import dqflags
from . import cube_build_wcs_util
from . import cube_overlap
#from . import cube_cloud_quick
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

        self.new_code = 0
        self.input_filenames = input_filenames
        self.pipeline = pipeline

        self.input_models = input_models # needed when building single mode IFU cubes
        self.output_name_base = output_name_base
        self.num_files = None

        self.instrument = instrument
        self.detector = detector
        self.list_par1 = list_par1
        self.list_par2 = list_par2

        self.instrument_info = instrument_info #dictionary class imported in cube_build.py
        self.master_table = master_table
        self.output_type = output_type
        self.scale1 = pars_cube.get('scale1')
        self.scale2 = pars_cube.get('scale2')
        self.scalew = pars_cube.get('scalew')
        self.rois = pars_cube.get('rois')
        self.roiw = pars_cube.get('roiw')
        self.spatial_size = None
        self.spectral_size = None
        self.interpolation = pars_cube.get('interpolation')
        self.coord_system = pars_cube.get('coord_system')
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

        self.soft_rad = None
        self.linear_wavelength = True
        self.roiw_table = None
        self.rois_table = None
        self.softrad_table = None
        self.weight_power_table = None
        self.wavelength_table = None

        self.cdelt1 = None
        self.cdelt2 = None
        self.cdelt3 = None
        self.crpix1 = None
        self.crpix2 = None
        self.crpix3 = None
        self.crval1 = None
        self.crval2 = None
        self.crval3 = None
        self.naxis1 = None
        self.naxis2 = None
        self.naxis3 = None
        self.cdelt3_normal = None

        self.a_min = 0
        self.a_max = 0
        self.b_min = 0
        self.b_max = 0
        self.lambda_min = 0
        self.lambda_max = 0
        self.xcoord = None
        self.ycoord = None
        self.zcoord = None

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
        major problems - cube_build is  being run incorrectly.

        """
        num1 = len(self.list_par1)
        num_files = 0
        for i in range(num1):
            this_a = self.list_par1[i]
            this_b = self.list_par2[i]
            n = len(self.master_table.FileMap[self.instrument][this_a][this_b])
            num_files = num_files + n
        self.num_files = num_files
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
#********************************************************************************



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
#        log.info('Ra average %f12.8', ra_ave)

        self.crval1 = ra_ave
        self.crval2 = dec_ave
        xi_center, eta_center = coord.radec2std(self.crval1, self.crval2,
                                                ra_ave, dec_ave)
        xi_min, eta_min = coord.radec2std(self.crval1, self.crval2, ra_min, dec_min)
        xi_max, eta_max = coord.radec2std(self.crval1, self.crval2, ra_max, dec_max)
#________________________________________________________________________________
# find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
# to find location of center abs of min values is how many pixels
        n1a = int(math.ceil(math.fabs(xi_min) / self.cdelt1))
        n2a = int(math.ceil(math.fabs(eta_min) / self.cdelt2))

        n1b = int(math.ceil(math.fabs(xi_max) / self.cdelt1))
        n2b = int(math.ceil(math.fabs(eta_max) / self.cdelt2))

        xi_min = 0.0 - (n1a * self.cdelt1) - (self.cdelt1 / 2.0)
        xi_max = (n1b * self.cdelt1) + (self.cdelt1 / 2.0)

        eta_min = 0.0 - (n2a * self.cdelt2) - (self.cdelt2 / 2.0)
        eta_max = (n2b * self.cdelt2) + (self.cdelt2 / 2.0)

        self.crpix1 = float(n1a) + 1.0
        self.crpix2 = float(n2a) + 1.0

        self.naxis1 = n1a + n1b
        self.naxis2 = n2a + n2b

        self.a_min = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max

# center of spaxels
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.cdelt1 / 2.0
        for i in range(self.naxis1):
            self.xcoord[i] = xstart
            xstart = xstart + self.cdelt1

        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.cdelt2 / 2.0
        for i in range(self.naxis2):
            self.ycoord[i] = ystart
            ystart = ystart + self.cdelt2

        ygrid = np.zeros(self.naxis2 * self.naxis1)
        xgrid = np.zeros(self.naxis2 * self.naxis1)

        k = 0
        ystart = self.ycoord[0]
        for i in range(self.naxis2):
            xstart = self.xcoord[0]
            for j in range(self.naxis1):
                xgrid[k] = xstart
                ygrid[k] = ystart
                xstart = xstart + self.cdelt1
                k = k + 1
            ystart = ystart + self.cdelt2

#        ycube,xcube = np.mgrid[0:self.naxis2,
#                               0:self.naxis1]
#        xcube = xcube.flatten()
#        ycube = ycube.flatten()

        self.xcenters = xgrid
        self.ycenters = ygrid
#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube
        self.cdelt3_normal = None
        if self.linear_wavelength:
            self.lambda_min = lambda_min
            self.lambda_max = lambda_max
            range_lambda = self.lambda_max - self.lambda_min
            self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))

         # adjust max based on integer value of naxis3
            lambda_center = (self.lambda_max + self.lambda_min) / 2.0
            self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.cdelt3
            self.lambda_max = self.lambda_min + (self.naxis3) * self.cdelt3

            self.zcoord = np.zeros(self.naxis3)
            self.crval3 = self.lambda_min
            self.crpix3 = 1.0
            zstart = self.lambda_min + self.cdelt3 / 2.0
            for i in range(self.naxis3):
                self.zcoord[i] = zstart
                zstart = zstart + self.cdelt3

        else:
            self.naxis3 = len(self.wavelength_table)
            self.zcoord = np.asarray(self.wavelength_table)
            self.crval3 = self.wavelength_table[0]
            self.crpix3 = 1.0
# set up the cdelt3_normal normalizing array used in cube_cloud.py
        cdelt3_normal = np.zeros(self.naxis3)
        for j in range(self.naxis3 - 1):
            cdelt3_normal[j] = self.zcoord[j + 1] - self.zcoord[j]

        cdelt3_normal[self.naxis3 - 1] = cdelt3_normal[self.naxis3 - 2]
        self.cdelt3_normal = cdelt3_normal
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
        self.naxis1 = int(math.ceil(range_a / self.cdelt1))

        # adjust min and max based on integer value of naxis1
        a_center = (self.a_max + self.a_min) / 2.0
        self.a_min = a_center - (self.naxis1 / 2.0) * self.cdelt1
        self.a_max = a_center + (self.naxis1 / 2.0) * self.cdelt1

        self.xcoord = np.zeros(self.naxis1)
        self.crval1 = self.a_min
        self.crpix1 = 0.5
        xstart = self.a_min + self.cdelt1 / 2.0
        for i in range(self.naxis1):
            self.xcoord[i] = xstart
            xstart = xstart + self.cdelt1
#_______________________________________________________________________
        #set up the lambda (z) coordinate of the cube
        range_lambda = self.lambda_max - self.lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))

         # adjust max based on integer value of naxis3
        lambda_center = (self.lambda_max + self.lambda_min) / 2.0

        self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.cdelt3
        self.lambda_max = lambda_center + (self.naxis3 / 2.0) * self.cdelt3

        self.lambda_max = self.lambda_min + (self.naxis3) * self.cdelt3

        self.zcoord = np.zeros(self.naxis3)
        self.crval3 = self.lambda_min
        self.crpix3 = 1.0
        zstart = self.lambda_min + self.cdelt3 / 2.0

        for i in range(self.naxis3):
            self.zcoord[i] = zstart
            zstart = zstart + self.cdelt3
#_______________________________________________________________________
        # set up the naxis2 parameters
        range_b = self.b_max - self.b_min

        self.naxis2 = int(math.ceil(range_b / self.cdelt2))
        b_center = (self.b_max + self.b_min) / 2.0
    # adjust min and max based on integer value of naxis2
        self.b_max = b_center + (self.naxis2 / 2.0) * self.cdelt2
        self.b_min = b_center - (self.naxis2 / 2.0) * self.cdelt2

        self.ycoord = np.zeros(self.naxis2)
        self.crval2 = self.b_min
        self.crpix2 = 0.5
        ystart = self.b_min + self.cdelt2 / 2.0
        for i in range(self.naxis2):
            self.ycoord[i] = ystart
            ystart = ystart + self.cdelt2

#_______________________________________________________________________
    def print_cube_geometry(self):

        """
        Print out the general properties of the size of the IFU Cube
        """

        log.info('Cube Geometry:')
        if self.coord_system == 'alpha-beta':
            log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (alpha,beta arc sec)')
        else:
            log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(arc sec)  MIN & Max (xi,eta arc sec)')
            log.info('Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                     self.naxis1, self.crpix1, self.crval1, self.cdelt1,
                     self.a_min, self.a_max)
            log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                     self.naxis2, self.crpix2, self.crval2, self.cdelt2,
                     self.b_min, self.b_max)
            if self.linear_wavelength:
                log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(microns)  MIN & Max (microns)')
                log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                         self.naxis3, self.crpix3, self.crval3, self.cdelt3,
                         self.lambda_min, self.lambda_max)

            if not self.linear_wavelength:
                log.info('Non-linear wavelength dimension, CDELT3 variable')
                log.info('axis#  Naxis  CRPIX    CRVAL     MIN & Max (microns)')
                log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f',
                         self.naxis3, self.crpix3, self.crval3,
                         self.wavelength_table[0], self.wavelength_table[self.naxis3 - 1])

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
        1. Loop over every band contained in the IFU cube and read in the data
        associated with the band
        2. map_detector_to_output_frame: Maps the detector data to the cube output coordinate system
        3. For each mapped detector pixel the ifu cube spaxel located in the region of
        interest. There are three different routines to do this step each of them use
        a slighly different weighting function in how to combine the detector fluxs that
        fall within a region of influence from the spaxel center
        a. cube_cloud:match_det2_cube_msm: This routine uses the modified
        shepard method to determing the weighting function, which weights the detector
        fluxes based on the distance between the detector center and spaxel center.
        b. cube_cloud:match_det2_cube_miripsf the weighting function based  width of the
        psf and lsf.
        c. cube_overlap.match_det2cube is only for single exposure, single band cubes and
        the ifucube in created in the detector plane. The weighting function is based on
        the overlap of between the detector pixel and spaxel. This method is simplified
        to determine the overlap in the alpha-wavelength plane.
        4. find_spaxel_flux: find the final flux assoicated with each spaxel
        5. setup_ifucube
        6. output_ifucube

        Returns
        -------
        Returns an ifu cube

        """

        self.output_name = self.define_cubename()
        total_num = self.naxis1 * self.naxis2 * self.naxis3
        self.spaxel_flux = np.zeros(total_num)
        self.spaxel_weight = np.zeros(total_num)
        self.spaxel_iflux = np.zeros(total_num)

        spaxel_ra = None
        spaxel_dec = None
        spaxel_wave = None
#________________________________________________________________________________
# Only preformed if weighting = MIRIPSF, first convert xi,eta cube to
# v2,v3,wave. This information if past to cube_cloud and for each
# input_model the v2,v3, wave is converted to alpha,beta in detector plane
# ra,dec, wave is independent of input_model
# v2,v3, alpha,beta depends on the input_model
        if self.weighting == 'miripsf':
            spaxel_ra = np.zeros(total_num)
            spaxel_dec = np.zeros(total_num)
            spaxel_wave = np.zeros(total_num)

            nxy = self.xcenters.size
            nz = self.zcoord.size
            for iz in range(nz):
                istart = iz * nxy
                for ixy in range(nxy):
                    ii = istart + ixy
                    spaxel_ra[ii], spaxel_dec[ii] = coord.std2radec(self.crval1,
                                                                    self.crval2,
                                                                    self.xcenters[ixy],
                                                                    self.ycenters[ixy])
                    spaxel_wave[ii] = self.zcoord[iz]
#________________________________________________________________________________

        subtract_background = True

        # now need to loop over every file that covers this channel/subchannel (MIRI)
        # or Grating/filter(NIRSPEC)
        #and map the detector pixels to the cube spaxel

        number_bands = len(self.list_par1)
        t0 = time.time()
        for i in range(number_bands):
            this_par1 = self.list_par1[i]
            this_par2 = self.list_par2[i]
            nfiles = len(self.master_table.FileMap[self.instrument][this_par1][this_par2])
#________________________________________________________________________________
# loop over the files that cover the spectral range the cube is for
            for k in range(nfiles):
                ifile = self.master_table.FileMap[self.instrument][this_par1][this_par2][k]
                self.this_cube_filenames.append(ifile)
                log.debug("Working on Band defined by:%s %s ", this_par1, this_par2)
#--------------------------------------------------------------------------------
                if self.interpolation == 'pointcloud':
                    t0 = time.time()
                    pixelresult = self.map_detector_to_outputframe(this_par1,
                                                                   this_par2,
                                                                   subtract_background,
                                                                   ifile)

                    coord1, coord2, wave, flux, rois_pixel, roiw_pixel, weight_pixel,\
                        softrad_pixel, alpha_det, beta_det = pixelresult
                    t1 = time.time()
                    log.info("Time to transform pixels to output frame = %.1f.s" % (t1 - t0,))
                    if self.weighting == 'msm':
                        t0 = time.time()
#                        if self.new_code:
#                            print('calling new cube_cloud')
#                            cube_cloud_quick.match_det2cube_msm(self.naxis1,self.naxis2,self.naxis3,
#                                                              self.cdelt1,self.cdelt2,self.cdelt3,
#                                                              self.rois,self.roiw,self.weight_power,
#                                                              self.xcoord,self.ycoord,self.zcoord,
#                                                              self.spaxel_flux,
#                                                              self.spaxel_weight,
#                                                              self.spaxel_iflux,
#                                                              flux,
#                                                              coord1,coord2,wave)
#                        else:
#                            print('calling old cube_cloud')
                        cube_cloud.match_det2cube_msm(self.naxis1, self.naxis2, self.naxis3,
                                                      self.cdelt1, self.cdelt2,
                                                      self.cdelt3_normal,
                                                      self.xcenters, self.ycenters, self.zcoord,
                                                      self.spaxel_flux,
                                                      self.spaxel_weight,
                                                      self.spaxel_iflux,
                                                      flux,
                                                      coord1, coord2, wave,
                                                      rois_pixel, roiw_pixel, weight_pixel,
                                                      softrad_pixel)


                        t1 = time.time()
                        log.info("Time to match file to ifucube = %.1f.s" % (t1 - t0,))
#________________________________________________________________________________
                    elif self.weighting == 'miripsf':
                        with datamodels.IFUImageModel(ifile) as input_model:
                            wave_resol = self.instrument_info.Get_RP_ave_Wave(this_par1,
                                                                              this_par2)

                            alpha_resol = self.instrument_info.Get_psf_alpha_parameters()
                            beta_resol = self.instrument_info.Get_psf_beta_parameters()

                            worldtov23 = input_model.meta.wcs.get_transform("world", "v2v3")
                            v2ab_transform = input_model.meta.wcs.get_transform('v2v3',
                                                                'alpha_beta')

                            spaxel_v2, spaxel_v3, zl = worldtov23(spaxel_ra,
                                                                  spaxel_dec,
                                                                  spaxel_wave)

                            spaxel_alpha, spaxel_beta, spaxel_wave = v2ab_transform(spaxel_v2,
                                                                                    spaxel_v3,
                                                                                    zl)
                            cube_cloud.match_det2cube_miripsf(alpha_resol,
                                                              beta_resol,
                                                              wave_resol,
                                                              self.naxis1, self.naxis2, self.naxis3,
                                                              self.xcenters, self.ycenters, self.zcoord,
                                                              self.spaxel_flux,
                                                              self.spaxel_weight,
                                                              self.spaxel_iflux,
                                                              spaxel_alpha, spaxel_beta, spaxel_wave,
                                                              flux,
                                                              coord1, coord2, wave,
                                                              alpha_det, beta_det,
                                                              rois_pixel, roiw_pixel, weight_pixel,
                                                              softrad_pixel)
#--------------------------------------------------------------------------------
#2D area method - only works for single files and coord_system = 'alpha-beta'
#--------------------------------------------------------------------------------
                elif self.interpolation == 'area':
                    with datamodels.IFUImageModel(ifile) as input_model:
                        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                              'alpha_beta')
                        start_region = self.instrument_info.GetStartSlice(this_par1)
                        end_region = self.instrument_info.GetEndSlice(this_par1)
                        regions = list(range(start_region, end_region + 1))
                        t0 = time.time()
                        for i in regions:
                            log.info('Working on Slice # %d', i)
                            y, x = (det2ab_transform.label_mapper.mapper == i).nonzero()

# getting pixel corner - ytop = y + 1 (routine fails for y = 1024)
                            index = np.where(y < 1023)
                            y = y[index]
                            x = x[index]
                            cube_overlap.match_det2cube(x, y, i,
                                                        start_region,
                                                        input_model,
                                                        det2ab_transform,
                                                        self.spaxel_flux,
                                                        self.spaxel_weight,
                                                        self.spaxel_iflux,
                                                        self.xcoord, self.zcoord,
                                                        self.crval1, self.crval3,
                                                        self.cdelt1, self.cdelt3,
                                                        self.naxis1, self.naxis2)
                        t1 = time.time()

                        log.info("Time Map All slices on Detector to Cube = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
# Mapped all data to cube or Point Cloud
# now determine Cube Spaxel flux

        t0 = time.time()
        self.find_spaxel_flux()

        t1 = time.time()
        log.info("Time to find Cube Flux= %.1f.s" % (t1 - t0,))

        ifucube_model = self.setup_ifucube(0)
#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube
        self.update_ifucube(ifucube_model)
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

        for j in range(n):
            log.info("Working on next Single IFU Cube  = %i" % (j + 1))
            t0 = time.time()
# for each new data model create a new spaxel
            total_num = self.naxis1 * self.naxis2 * self.naxis3
            self.spaxel_flux = np.zeros(total_num)
            self.spaxel_weight = np.zeros(total_num)
            self.spaxel_iflux = np.zeros(total_num)

            subtract_background = False

            pixelresult = self.map_detector_to_outputframe(this_par1, this_par2,
                                                           subtract_background,
                                                           self.input_models[j])

            coord1, coord2, wave, flux, rois_pixel, roiw_pixel, weight_pixel, \
                softrad_pixel, alpha_det, beta_det = pixelresult

            cube_cloud.match_det2cube_msm(self.naxis1, self.naxis2, self.naxis3,
                                          self.cdelt1, self.cdelt2,
                                          self.cdelt3_normal,
                                          self.xcenters, self.ycenters, self.zcoord,
                                          self.spaxel_flux,
                                          self.spaxel_weight,
                                          self.spaxel_iflux,
                                          flux,
                                          coord1, coord2, wave,
                                          rois_pixel, roiw_pixel, weight_pixel,
                                          softrad_pixel)
#_______________________________________________________________________
# shove Flux and iflux in the  final ifucube
            self.find_spaxel_flux()
# now determine Cube Spaxel flux

            ifucube_model = self.setup_ifucube(j)
            self.update_ifucube(ifucube_model)
            t1 = time.time()
            log.info("Time Create Single ifucube  = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
            single_ifucube_container.append(ifucube_model)

        return single_ifucube_container
#********************************************************************************

    def determine_cube_parameters(self):
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
        # initialize
        wave_roi = None
        weight_power = None

        number_bands = len(self.list_par1)
        spaxelsize = np.zeros(number_bands)
        spectralsize = np.zeros(number_bands)
        rois = np.zeros(number_bands)
        roiw = np.zeros(number_bands)
        power = np.zeros(number_bands)
        softrad = np.zeros(number_bands)
        minwave = np.zeros(number_bands)
        maxwave = np.zeros(number_bands)

        for i in range(number_bands):
            if self.instrument == 'MIRI':
                par1 = self.list_par1[i]
                par2 = self.list_par2[i]
            elif self.instrument == 'NIRSPEC':
                par1 = self.list_par1[i]
                par2 = self.list_par2[i]

            roiw[i] = self.instrument_info.GetWaveRoi(par1, par2)
            rois[i] = self.instrument_info.GetSpatialRoi(par1, par2)

            a_scale, b_scale, w_scale = self.instrument_info.GetScale(par1,
                                                                      par2)
            spaxelsize[i] = a_scale
            spectralsize[i] = w_scale

            power[i] = self.instrument_info.GetMSMPower(par1, par2)
            softrad[i] = self.instrument_info.GetSoftRad(par1, par2)
            minwave[i] = self.instrument_info.GetWaveMin(par1, par2)
            maxwave[i] = self.instrument_info.GetWaveMax(par1, par2)
# Check the spatial size. If it is the same for the array set up the parameters
        all_same = np.all(spaxelsize == spaxelsize[0])
        if all_same:
            self.spatial_size = spaxelsize[0]
            spatial_roi = rois[0]
        else:
            index_min = np.argmin(spaxelsize)
            self.spatial_size = spaxelsize[index_min]
            spatial_roi = rois[index_min]
# find min and max wavelength
        min_wave = np.amin(minwave)
        max_wave = np.amax(maxwave)

        if self.wavemin is None:
            self.wavemin = min_wave
        else:
            self.wavemin = np.float64(self.wavemin)

        if self.wavemax is None:
            self.wavemax = max_wave
        else:
            self.wavemax = np.float64(self.wavemax)

# now check spectral step
        all_same_spectral = np.all(spectralsize == spectralsize[0])
        if all_same_spectral:
            self.spectral_size = spectralsize[0]
            wave_roi = roiw[0]
            self.soft_rad = softrad[0]
            weight_power = power[0]
        else:
            self.linear_wavelength = False
            if self.instrument == 'MIRI':
                table = self.instrument_info.Get_multichannel_table()
                table_wavelength, table_sroi, table_wroi, table_power, table_softrad = table

# getting MIRI Table Values
            elif self.instrument == 'NIRSPEC':
# determine if have Prism, Medium or High resolution
                med = ['g140m', 'g235m', 'g395m']
                high = ['g140h', 'g235h', 'g395h']
                prism = ['prism']

                for i in range(number_bands):
                    par1 = self.list_par1[i]
                    if par1 in prism:
                        table = self.instrument_info.Get_prism_table()
                    if par1 in med:
                        table = self.instrument_info.Get_med_table()
                    if par1 in high:
                        table = self.instrument_info.Get_high_table()
                    table_wavelength, table_sroi, table_wroi, table_power, table_softrad = table
# based on Min and Max wavelength - pull out the tables values that fall in this range
            # find the closest table entries to the self.wavemin and self.wavemax limits
            imin = (np.abs(table_wavelength - self.wavemin)).argmin()
            imax = (np.abs(table_wavelength - self.wavemax)).argmin()
            if imin > 1 and table_wavelength[imin] > self.wavemin: imin = imin - 1
            if (imax < len(table_wavelength) and
                self.wavemax > table_wavelength[imax]): imax = imax + 1
#            print('index of wavelength values',imin,imax)

            self.roiw_table = table_wroi[imin:imax]
            self.rois_table = table_sroi[imin:imax]
            self.softrad_table = table_softrad[imin:imax]
            self.weight_power_table = table_power[imin:imax]
            self.wavelength_table = table_wavelength[imin:imax]

# check if the user has set the cube parameters to use
        if self.rois == 0: self.rois = spatial_roi
        if self.output_type == 'single' or self.num_files < 4:
            self.rois = self.rois * 1.5
            log.info('Increasing spatial region of interest' +
                     'default value set for 4 dithers %f', self.rois)
        if self.scale1 != 0:
            self.spatial_size = self.scale1
        if self.scalew != 0:
            self.spectral_size = self.scalew
            self.linear_wavelength = True
            #set wave_roi, weight_power, soft_rad to same values if they are in  list
        if self.roiw == 0: self.roiw = wave_roi
        if self.weight_power == 0: self.weight_power = weight_power

#        print('spatial size', self.spatial_size)
#        print('spectral size', self.spectral_size)
#        print('spatial roi', self.rois)
#        print('wave min and max', self.wavemin, self.wavemax)
#        print('linear wavelength', self.linear_wavelength)
#        print('roiw', self.roiw)



#        if self.interpolation == 'pointcloud':
#            log.info('Region of interest spatial, wavelength  %f %f', self.rois, self.roiw)

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

        If the coordinate system is alpha-beta (MIRI) then min and max
        coordinates of alpha (arc sec), beta (arc sec) and lambda (microns)
        If the coordinate system is world then the min and max of
        ra(degress), dec (degrees) and lambda (microns) is returned.

        """
#________________________________________________________________________________
        self.cdelt1 = self.spatial_size
        self.cdelt2 = self.spatial_size
        if self.linear_wavelength:
            self.cdelt3 = self.spectral_size

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
            for k in range(n):
                amin = 0.0
                amax = 0.0
                bmin = 0.0
                bmax = 0.0
                lmin = 0.0
                lmax = 0.0

                ifile = self.master_table.FileMap[self.instrument][this_a][this_b][k]
#________________________________________________________________________________
# Open the input data model
# Find the footprint of the image
                with datamodels.IFUImageModel(ifile) as input_model:
                    if self.instrument == 'NIRSPEC':
                        flag_data = 0
#                        t0 = time.time()
                        ch_footprint = cube_build_wcs_util.find_footprint_NIRSPEC(
                            input_model,
                            flag_data,
                            self.coord_system)
#                        t1 = time.time()
#                        print('time to find footprint',t1-t0)

                        amin, amax, bmin, bmax, lmin, lmax = ch_footprint
#                        amin, bmin, amax, bmin2, amax2, bmax, amin2, bmax2 = \
#                            input_model.meta.wcsinfo.s_region
#                        print(amin, amax, bmin, bmax, amin2, amax2, bmin2, bmax2)
#________________________________________________________________________________
                    if self.instrument == 'MIRI':
                        ch_footprint = cube_build_wcs_util.find_footprint_MIRI(
                            input_model,
                            this_a,
                            self.instrument_info,
                            self.coord_system)
                        amin, amax, bmin, bmax, lmin, lmax = ch_footprint
#                        print(amin,amax,bmin,bmax)
#                        region_footprint = input_model.meta.wcsinfo.s_region
#                        print('region footprint', region_footprint)
#                        print(type(region_footprint))
#                        amin, bmin, amax, bmin2, amax2, bmax, amin2, bmax2 = region_footprint
#                        print(amin, amax, bmin, bmax, amin2, amax2, bmin2, bmax2)
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

        print('wavelengths', final_lambda_min, final_lambda_max, self.wavemin, self.wavemax)

        final_lambda_min = self.wavemin
        final_lambda_max = self.wavemax
#________________________________________________________________________________
        if self.instrument == 'MIRI' and self.coord_system == 'alpha-beta':
        #  we have a 1 to 1 mapping in beta dimension.
            nslice = self.instrument_info.GetNSlice(parameter1[0])
            log.info('Beta Scale %f ', self.cdelt2)
            self.cdelt2 = (final_b_max - final_b_min) / nslice
            final_b_max = final_b_min + (nslice) * self.cdelt2
            log.info('Changed the Beta Scale dimension so we have 1 -1 mapping between beta and slice #')
            log.info('New Beta Scale %f ', self.cdelt2)
#________________________________________________________________________________
# Test that we have data (NIRSPEC NRS2 only has IFU data for 3 configurations)
        test_a = final_a_max - final_a_min
        test_b = final_b_max - final_b_min
        tolerance1 = 0.00001
        if(test_a < tolerance1 or test_b < tolerance1):
            log.info('No Valid IFU slice data found %f %f ', test_a, test_b)
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
    def map_detector_to_outputframe(self, this_par1, this_par2,
                                    subtract_background,
                                    ifile):
        from ..mrs_imatch.mrs_imatch_step import apply_background_2d
#********************************************************************************
        """
        Short Summary
        -------------
        Loop over a file and map the detector pixels to the output frame of the ifcube

        Return the coordinates of the pixel in the output frame as well as the associatied
        flux for those pixels.

        Parameter
        ----------

        Returns
        -------
        if(interpolation = area - only valid for alpha-beta
        or
        if(interpolation = pointcloud
        """

# intitalize alpha_det and beta_det to None. These are filled in if the instrument
# is MIRI and the weighting is miripsf

        alpha_det = None
        beta_det = None
        coord1 = None
        coord2 = None
        flux = None
        wave = None
# Open the input data model
        with datamodels.IFUImageModel(ifile) as input_model:
            # check if background sky matching as been done
            # mrs_imatch step. THis is only for MRS data at this time
            # but go head and check it before splitting by instrument
            # the polynomial should be empty for NIRSPEC

            num_ch_bgk = len(input_model.meta.background.polynomial_info)
            if(num_ch_bgk > 0 and subtract_background):
                for ich_num in range(num_ch_bgk):
                    poly = input_model.meta.background.polynomial_info[ich_num]
                    poly_ch = poly.channel
                    if(poly_ch == this_par1):
                        apply_background_2d(input_model, poly_ch, subtract=True)

#********************************************************************************
#--------------------------------------------------------------------------------
            if self.instrument == 'MIRI':
                xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_par1)
                y, x = np.mgrid[:1024, xstart:xend]
                y = np.reshape(y, y.size)
                x = np.reshape(x, x.size)
                if self.coord_system == 'world':
                    ra, dec, wave = input_model.meta.wcs(x, y)
                    valid1 = ~np.isnan(ra)
                    ra = ra[valid1]
                    dec = dec[valid1]
                    wave = wave[valid1]
                    x = x[valid1]
                    y = y[valid1]
                    if self.weighting == 'miripsf':
                        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                          'alpha_beta')
                        alpha, beta, lam = det2ab_transform(x, y)
                elif self.coord_system == 'alpha-beta':
                    det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                          'alpha_beta')
                    alpha, beta, wave = det2ab_transform(x, y)
                    valid1 = ~np.isnan(coord1)
                    alpha = alpha[valid1]
                    beta = beta[valid1]
                    wave = wave[valid1]
                    x = x[valid1]
                    y = y[valid1]
#________________________________________________________________________________
            elif self.instrument == 'NIRSPEC':
                # initialize the ra,dec, and wavelength arrays
                # we will loop over slices and fill in values
                # the flag_det will be set when a slice pixel is filled in
                #   at the end we will use this flag to pull out valid data
                ra_det = np.zeros((2048, 2048))
                dec_det = np.zeros((2048, 2048))
                lam_det = np.zeros((2048, 2048))
                flag_det = np.zeros((2048, 2048))
                # for NIRSPEC each file has 30 slices
                # wcs information access seperately for each slice
                nslices = 30
                log.info("Mapping each NIRSPEC slice to sky, this takes a while for NIRSPEC data")
                for ii in range(nslices):
                    slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)
                    x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
                    ra, dec, lam = slice_wcs(x, y)

                    # the slices are curved on detector so a rectangular region
                    # returns NaNs
                    valid = ~np.isnan(lam)
                    ra = ra[valid]
                    dec = dec[valid]
                    lam = lam[valid]
                    x = x[valid]
                    y = y[valid]

                    xind = _toindex(x)
                    yind = _toindex(y)
                    xind = np.ndarray.flatten(xind)
                    yind = np.ndarray.flatten(yind)
                    ra = np.ndarray.flatten(ra)
                    dec = np.ndarray.flatten(dec)
                    lam = np.ndarray.flatten(lam)
                    ra_det[yind, xind] = ra
                    dec_det[yind, xind] = dec
                    lam_det[yind, xind] = lam
                    flag_det[yind, xind] = 1
                    # done looping  - pull out valid values
                valid_data = np.where(flag_det == 1)
                y, x = valid_data
                ra = ra_det[valid_data]
                dec = dec_det[valid_data]
                wave = lam_det[valid_data]
# ________________________________________________________________________________
# The following is for both MIRI and NIRSPEC
# grab the flux and DQ values for these pixles
            flux_all = input_model.data[y, x]
            dq_all = input_model.dq[y, x]
            valid2 = np.isfinite(flux_all)
#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded
# from the cube mapping
            all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] +
                         dqflags.pixel['NON_SCIENCE'] +
                         dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] +
                         dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])
            # find the location of all the values to reject in cube building
            good_data = np.where((np.bitwise_and(dq_all, all_flags) == 0) & (valid2 == True))
                    # good data holds the location of pixels we want to map to cube
            flux = flux_all[good_data]
            wave = wave[good_data]
# based on the wavelength define the sroi, wroi, weight_power and softrad to use
# in matching detector to spaxel values

            rois_det = np.zeros(wave.shape)
            roiw_det = np.zeros(wave.shape)
            weight_det = np.zeros(wave.shape)
            softrad_det = np.zeros(wave.shape)
            if self.linear_wavelength:
                rois_det[:] = self.rois
                roiw_det[:] = self.roiw
                weight_det[:] = self.weight_power
                softrad_det[:] = self.soft_rad
            else:
                #for each wavelength find the closest point in the self.wavelength_table
                for iw, w in enumerate(wave):
                    ifound = (np.abs(self.wavelength_table - w)).argmin()
                    rois_det[iw] = self.rois_table[ifound]
                    roiw_det[iw] = self.roiw_table[ifound]
                    softrad_det[iw] = self.softrad_table[ifound]
                    weight_det[iw] = self.weight_power_table[ifound]

            if self.coord_system == 'world':
                ra_use = ra[good_data]
                dec_use = dec[good_data]
                coord1, coord2 = coord.radec2std(self.crval1, self.crval2, ra_use, dec_use)
                if self.weighting == 'miripsf':
                    alpha_det = alpha[good_data]
                    beta_det = beta[good_data]
            elif self.coord_system == 'alpha-beta':
                coord1 = alpha[good_data]
                coord2 = beta[good_data]

        return coord1, coord2, wave, flux, rois_det, roiw_det, weight_det, \
            softrad_det, alpha_det, beta_det

#********************************************************************************
    def find_spaxel_flux(self):
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
# currently these are the same but in the future there could be a difference in
# how the spaxel flux is determined according to self.interpolation.

        if self.interpolation == 'area':
            good = self.spaxel_iflux > 0
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]
        elif self.interpolation == 'pointcloud':
            good = self.spaxel_iflux > 0
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]

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

        if self.linear_wavelength:
            ifucube_model = datamodels.IFUCubeModel(data=data, dq=dq_cube,
                                                    err=err_cube,
                                                    weightmap=idata)
        else:
            nelem = np.array([])
            wave = np.asarray(self.wavelength_table, dtype=np.float32)
            # for now we need to pad wavelength to fix datamodel size
            num = len(wave)
            nn = 2000 # This number needs to match the shape of the
            # wavetable['wavelength'] value in the ifucube.schema.
            # In the future it would be good to remove that we have to
            # set the size of this value in the schema.

            wave = np.pad(wave, (0, nn - num), 'constant')
            newnum = len(wave)
            allwave = np.zeros((1, 1, newnum))
            allwave[0, 0, :] = wave

            nelem = np.append(nelem, num)
# to get the data in the correct format (an array in a single cell in the fit table)
# I had to zip data. 
            alldata = np.array(list(zip(np.array(nelem), np.array(allwave))),
                               dtype=datamodels.IFUCubeModel().wavetable.dtype)


            ifucube_model = datamodels.IFUCubeModel(data=data, dq=dq_cube,
                                                    err=err_cube,
                                                    weightmap=idata,
                                                    wavetable=alldata)



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
# fill in Channel for MIRI
        if self.instrument == 'MIRI':
            # fill in Channel output meta data
            num_ch = len(self.list_par1)
            outchannel = self.list_par1[0]
            for m in range(1, num_ch):
                outchannel = outchannel + str(self.list_par1[m])

            outchannel = "".join(set(outchannel))
            outchannel = "".join(sorted(outchannel))
            ifucube_model.meta.instrument.channel = outchannel
            log.info('IFUChannel %s', ifucube_model.meta.instrument.channel)

#______________________________________________________________________
        ifucube_model.meta.wcsinfo.crval1 = self.crval1
        ifucube_model.meta.wcsinfo.crval2 = self.crval2
        ifucube_model.meta.wcsinfo.crpix1 = self.crpix1
        ifucube_model.meta.wcsinfo.crpix2 = self.crpix2

        ifucube_model.meta.wcsinfo.cdelt1 = self.cdelt1 / 3600.0
        ifucube_model.meta.wcsinfo.cdelt2 = self.cdelt2 / 3600.0
        if self.linear_wavelength:
            ifucube_model.meta.wcsinfo.crval3 = self.crval3
            ifucube_model.meta.wcsinfo.cdelt3 = self.cdelt3
            ifucube_model.meta.wcsinfo.ctype3 = 'WAVE'
            ifucube_model.meta.wcsinfo.crpix3 = self.crpix3
        else:
            ifucube_model.meta.wcsinfo.ctype3 = 'WAVE-TAB'
            ifucube_model.meta.wcsinfo.ps3_0 = 'WSC-TABLE'
            ifucube_model.meta.wcsinfo.ps3_1 = 'wavelength'

        ifucube_model.meta.wcsinfo.ctype1 = 'RA---TAN'
        ifucube_model.meta.wcsinfo.ctype2 = 'DEC--TAN'
        ifucube_model.meta.wcsinfo.cunit1 = 'deg'
        ifucube_model.meta.wcsinfo.cunit2 = 'deg'

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
        # weight_power is needed for single cubes. Linear Wavelengths
        # if non-linear wavelengths then this will be None
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
    def update_ifucube(self, ifucube_model):
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
    #pull out data into array and assign to ifucube data model
        temp_flux = self.spaxel_flux.reshape((self.naxis3, self.naxis2, self.naxis1))
        temp_wmap = self.spaxel_iflux.reshape((self.naxis3, self.naxis2, self.naxis1))

        ifucube_model.data = temp_flux
        ifucube_model.weightmap = temp_wmap
        ifucube_model.meta.cal_step.cube_build = 'COMPLETE'

#********************************************************************************
    def blend_output_metadata(self, IFUCube):

        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = IFUCube.meta.filename
        blendmeta.blendmodels(IFUCube, inputs=self.input_models,
                              output=output_file)
class IncorrectInput(Exception):
    pass

class AreaInterpolation(Exception):
    pass
