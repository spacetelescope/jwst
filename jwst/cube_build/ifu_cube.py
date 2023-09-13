""" Work horse routines used for building ifu spectra cubes
(including the main loop over files and the construction of
final spaxel fluxes)
"""

import numpy as np
import logging
import math

from astropy.stats import circmean
from astropy import units as u
from gwcs import wcstools

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stdatamodels.jwst.transforms.models import _toindex

from ..model_blender import blendmeta
from ..assign_wcs import pointing
from jwst.datamodels import ModelContainer
from ..assign_wcs import nirspec
from ..assign_wcs.util import wrap_ra
from . import cube_build_wcs_util
from . import cube_internal_cal
from . import coord
from ..mrs_imatch.mrs_imatch_step import apply_background_2d
from .cube_match_sky_pointcloud import cube_wrapper  # c extension
from .cube_match_sky_driz import cube_wrapper_driz  # c extension

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class IFUCubeData():

    def __init__(self,
                 pipeline,
                 input_models,
                 output_name_base,
                 output_type,
                 instrument,
                 list_par1,
                 list_par2,
                 instrument_info,
                 master_table,
                 **pars_cube):
        """ Class IFUCube holds the high level data for each IFU Cube
        """
        self.input_models_this_cube = []  # list of files use to make cube working on

        self.pipeline = pipeline

        self.input_models = input_models  # needed when building single mode IFU cubes
        self.output_name_base = output_name_base
        self.num_files = None

        self.instrument = instrument
        self.list_par1 = list_par1
        self.list_par2 = list_par2

        self.instrument_info = instrument_info  # dictionary class imported in cube_build.py
        self.master_table = master_table
        self.output_type = output_type

        self.scalexy = pars_cube.get('scalexy')
        self.scalew = pars_cube.get('scalew')
        self.ra_center = pars_cube.get('ra_center')
        self.dec_center = pars_cube.get('dec_center')
        self.cube_pa = pars_cube.get('cube_pa')
        self.nspax_x = pars_cube.get('nspax_x')
        self.nspax_y = pars_cube.get('nspax_y')
        self.rois = pars_cube.get('rois')
        self.roiw = pars_cube.get('roiw')
        self.debug_spaxel = pars_cube.get('debug_spaxel')

        self.spaxel_x, self.spaxel_y, self.spaxel_z = [int(val) for val in self.debug_spaxel.split()]
        self.spatial_size = None
        self.spectral_size = None
        self.interpolation = pars_cube.get('interpolation')
        self.coord_system = pars_cube.get('coord_system')
        self.wavemin = pars_cube.get('wavemin')
        self.wavemax = pars_cube.get('wavemax')
        self.weighting = pars_cube.get('weighting')
        self.weight_power = pars_cube.get('weight_power')
        self.skip_dqflagging = pars_cube.get('skip_dqflagging')
        self.suffix = pars_cube.get('suffix')
        self.num_bands = 0
        self.output_name = ''

        self.wavemin_user = False  # Check for NIRSpec if user has set wavelength limits
        self.wavemax_user = False
        self.soft_rad = None
        self.scalerad = None
        self.linear_wavelength = True
        self.roiw_table = None
        self.rois_table = None
        self.softrad_table = None
        self.scalerad_table = None
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
        self.rot_angle = None  # rotation angle between Ra-Dec and IFU local instrument plane

        self.a_min = 0
        self.a_max = 0
        self.b_min = 0
        self.b_max = 0
        self.lambda_min = 0
        self.lambda_max = 0
        self.xcoord = None
        self.ycoord = None
        self.zcoord = None

        self.tolerance_dq_overlap = 0.05  # spaxel has to have 5% overlap to flag in FOV
        self.overlap_partial = 4  # intermediate flag
        self.overlap_full = 2    # intermediate flag
        self.overlap_hole = dqflags.pixel['DO_NOT_USE']
        self.overlap_no_coverage = dqflags.pixel['NON_SCIENCE']

    # **************************************************************
    def check_ifucube(self):
        """ Perform some quick checks that the type of cube to be produced
        conforms to rules

        Raises
        ------
        IncorrectInput
          Interpolation = area was selected for when input data is more than
          one file or model

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
        if self.coord_system == "internal_cal":
            if num_files > 1:
                raise IncorrectInput("Cubes built in internal_cal coordinate system" +
                                     " are built from a single file, not multiple exposures")
            if len(self.list_par1) > 1:
                raise IncorrectInput("Only a single channel or grating " +
                                     " can be used to create cubes in internal_cal coordinate system." +
                                     " Use --output_type=band")

    # ________________________________________________________________________________
    def define_cubename(self):
        """ Define the base output name
        """
        if self.pipeline == 2:
            newname = self.output_name_base + self.suffix + '.fits'
        else:
            if self.instrument == 'MIRI':

                # Check to see if the output base name already contains the
                # field "clear", which sometimes shows up in IFU product
                # names created by the ASN rules. If so, strip it off, so
                # that the remaining suffixes created below form the entire
                # list of optical elements in the final output name.
                suffix = self.output_name_base[self.output_name_base.rfind('_') + 1:]
                if suffix in ['clear']:
                    self.output_name_base = self.output_name_base[:self.output_name_base.rfind('_')]

                # Now compose the appropriate list of optical element suffix names
                # based on MRS channel and sub-channel
                channels = []
                for ch in self.list_par1:
                    if ch not in channels:
                        channels.append(ch)
                    number_channels = len(channels)
                    ch_name = '_ch'
                    for i in range(number_channels):
                        ch_name = ch_name + channels[i]
                        if i < number_channels - 1:
                            ch_name = ch_name + '-'

                # Sort by inverse alphabetical, e.g. short -> medium -> long
                subchannels = sorted(list(set(self.list_par2)))[::-1]
                log.info(f"Subchannel listing: {subchannels}")
                number_subchannels = len(subchannels)
                b_name = ''
                for i in range(number_subchannels):
                    b_name = b_name + subchannels[i]
                b_name = b_name.lower()
                newname = self.output_name_base + ch_name + '-' + b_name
                if self.coord_system == 'internal_cal':
                    newname = self.output_name_base + ch_name + '-' + b_name + '_internal'
                if self.output_type == 'single':
                    newname = self.output_name_base + ch_name + '-' + b_name + '_single'
            # ________________________________________________________________________________
            elif self.instrument == 'NIRSPEC':

                # Check to see if the output base name already has a grating/prism
                # suffix attached. If so, strip it off, and let the following logic
                # add all necessary grating and filter suffixes.
                suffix = self.output_name_base[self.output_name_base.rfind('_') + 1:]
                if suffix in ['g140m', 'g235m', 'g395m', 'g140h', 'g235h', 'g395h', 'prism']:
                    self.output_name_base = self.output_name_base[:self.output_name_base.rfind('_')]

                fg_name = '_'
                for i in range(len(self.list_par1)):
                    fg_name = fg_name + self.list_par1[i] + '-' + self.list_par2[i]
                    if i < self.num_bands - 1:
                        fg_name = fg_name + '-'
                fg_name = fg_name.lower()
                newname = self.output_name_base + fg_name
                if self.output_type == 'single':
                    newname = self.output_name_base + fg_name + '_single'
                if self.coord_system == 'internal_cal':
                    newname = self.output_name_base + fg_name + '_internal'
        # ______________________________________________________________________________
        if self.output_type != 'single':
            log.info(f'Output Name: {newname}')
        return newname
    # _______________________________________________________________________

    def set_geometry(self, corner_a, corner_b, lambda_min, lambda_max):
        """ Based on the ra,dec and wavelength footprint set up the size
        of the cube in the tangent plane projected coordinate system.


        Parameters
        ----------
        footprint: tuple
          holds min and max or ra,dec, and wavelength for the cube
          footprint
        """

        ra_min = np.min(corner_a)
        ra_max = np.max(corner_a)

        dec_min = np.min(corner_b)
        dec_max = np.max(corner_b)
        dec_ave = (dec_min + dec_max) / 2.0

        # we can not average ra values because of the convergence
        # of hour angles.
        ravalues = np.zeros(2)
        ravalues[0] = ra_min
        ravalues[1] = ra_max

        # astropy circmean assumes angles are in radians,
        # we have angles in degrees
        ra_ave = circmean(ravalues * u.deg).value

        if self.ra_center is not None:
            self.crval1 = self.ra_center
        else:
            self.crval1 = ra_ave

        if self.dec_center is not None:
            self.crval2 = self.dec_center
        else:
            self.crval2 = dec_ave

        rot_angle = self.rot_angle

        # find the 4 corners
        xi_corner = []
        eta_corner = []
        num = len(corner_a)

        for i in range(num):
            xi, eta = coord.radec2std(self.crval1, self.crval2,
                                      corner_a[i], corner_b[i], rot_angle)
            xi_corner.append(xi)
            eta_corner.append(eta)
        xi_min = min(xi_corner)
        xi_max = max(xi_corner)
        eta_min = min(eta_corner)
        eta_max = max(eta_corner)
        # ________________________________________________________________________________
        # find the CRPIX1 CRPIX2 - xi and eta centered at 0,0
        # to find location of center abs of min values is how many pixels
        # we want a symmetric cube centered on xi,eta = 0
        xilimit = max(np.abs(xi_min), np.abs(xi_max))
        etalimit = max(np.abs(eta_min), np.abs(eta_max))

        na = math.ceil(xilimit / self.cdelt1) + 1
        nb = math.ceil(etalimit / self.cdelt2) + 1

        # if the user set the nspax_x or nspax_y then redefine na, nb
        # it is assumed that both values are ODD numbers
        # We want the central pixel to be the tangent point with na/nb pixels on either
        # side of central pixel. 
        if self.nspax_x is not None:
            na = int(self.nspax_x/2)
        if self.nspax_y is not None:
            nb = int(self.nspax_y/2)

        xi_min = 0.0 - (na * self.cdelt1) - (self.cdelt1 / 2.0)
        xi_max = (na * self.cdelt1) + (self.cdelt1 / 2.0)

        eta_min = 0.0 - (nb * self.cdelt2) - (self.cdelt2 / 2.0)
        eta_max = (nb * self.cdelt2) + (self.cdelt2 / 2.0)

        self.crpix1 = float(na) + 1.0
        self.crpix2 = float(nb) + 1.0

        self.naxis1 = na * 2 + 1
        self.naxis2 = nb * 2 + 1

        self.a_min = xi_min
        self.a_max = xi_max
        self.b_min = eta_min
        self.b_max = eta_max
        # center of spaxels
        self.xcoord = np.zeros(self.naxis1)
        xstart = xi_min + self.cdelt1 / 2.0
        self.xcoord = np.arange(start=xstart, stop=xstart + self.naxis1 * self.cdelt1, step=self.cdelt1)

        self.ycoord = np.zeros(self.naxis2)
        ystart = eta_min + self.cdelt2 / 2.0
        self.ycoord = np.arange(start=ystart, stop=ystart + self.naxis2 * self.cdelt2, step=self.cdelt2)
        # depending on the naxis and cdelt values the x,ycoord can have 1 more element than naxis.
        # Clean up arrays dropping extra values at the end.
        self.xcoord = self.xcoord[0:self.naxis1]
        self.ycoord = self.ycoord[0:self.naxis2]

        xv, yv = np.meshgrid(self.xcoord, self.ycoord)
        self.xcenters = xv.flatten()
        self.ycenters = yv.flatten()
        # _______________________________________________________________________
        # set up the lambda (z) coordinate of the cube
        self.cdelt3_normal = None
        if self.linear_wavelength:
            self.lambda_min = lambda_min
            self.lambda_max = lambda_max
            range_lambda = self.lambda_max - self.lambda_min
            self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))

            # adjust max based on integer value of naxis3
            self.lambda_max = self.lambda_min + (self.naxis3) * self.cdelt3

            self.zcoord = np.zeros(self.naxis3)
            # CRPIX3 for FITS is 1 (center of first pixel)
            # CRVAL3 then is lambda_min + self.cdelt3/ 2.0, which is also zcoord[0]
            # Note that these are all values at the center of a spaxel
            self.crval3 = self.lambda_min + self.cdelt3 / 2.0
            self.crpix3 = 1.0
            zstart = self.lambda_min + self.cdelt3 / 2.0
            self.zcoord = np.arange(start=zstart, stop=self.lambda_max, step=self.cdelt3)
            self.zcoord = self.zcoord[0:self.naxis3]
        else:
            self.naxis3 = len(self.wavelength_table)
            self.zcoord = np.asarray(self.wavelength_table)
            self.crval3 = self.wavelength_table[0]
            self.crpix3 = 1.0
        # set up the cdelt3_normal normalizing array used
        cdelt3_normal = np.zeros(self.naxis3)
        for j in range(self.naxis3 - 1):
            cdelt3_normal[j] = self.zcoord[j + 1] - self.zcoord[j]

        cdelt3_normal[self.naxis3 - 1] = cdelt3_normal[self.naxis3 - 2]
        self.cdelt3_normal = cdelt3_normal
    # _______________________________________________________________________

    def set_geometryAB(self, corner_a, corner_b, lambda_min, lambda_max):
        """Based on the along slice, across slice and wavelength footprint set up the
        size of the cube in internal IFU plane.

        This will be a single exposure cube - small FOV assume
        rectangular coord system.

        Parameters
        ----------
        footprint : tuple
           Holds the min and max alpha, beta and wavelength values of
           cube on sky
        """

        self.a_min = np.min(corner_a)
        self.a_max = np.max(corner_a)

        self.b_min = np.min(corner_b)
        self.b_max = np.max(corner_b)
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max

        # along slice: a
        # across slice: b
        alimit = max(np.abs(self.a_min), np.abs(self.a_max))

        range_b = self.b_max - self.b_min
        if self.instrument == 'MIRI':
            # self.cdelt1 = self.cdelt2 # make cubes same scaling. MIRI EC team requested this removed (2/16/21)
            along_cdelt = self.cdelt1

            n1a = math.ceil(alimit / along_cdelt)
            n1b = math.ceil(alimit / along_cdelt)

            self.a_min = 0.0 - (n1a * along_cdelt) - (along_cdelt / 2.0)
            self.a_max = (n1b * along_cdelt) + (along_cdelt / 2.0)
            self.naxis1 = n1a + n1b + 1
            along_naxis = self.naxis1

            self.naxis2 = int(math.ceil(range_b / self.cdelt2))
            across_naxis = self.naxis2
            across_cdelt = self.cdelt2

        if self.instrument == 'NIRSPEC':
            along_cdelt = self.cdelt2

            n1a = math.ceil(alimit / along_cdelt)
            n1b = math.ceil(alimit / along_cdelt)

            self.a_min = 0.0 - (n1a * along_cdelt) - (along_cdelt / 2.0)
            self.a_max = (n1b * along_cdelt) + (along_cdelt / 2.0)

            self.naxis2 = n1a + n1b + 1
            along_naxis = self.naxis2
            self.naxis1 = int(math.ceil(range_b / self.cdelt1))
            across_naxis = self.naxis1
            across_cdelt = self.cdelt1

        acoord = np.zeros(along_naxis)
        astart = self.a_min + (along_cdelt / 2.0)
        for i in range(along_naxis):
            acoord[i] = astart
            astart = astart + along_cdelt

        # set up the across slice  parameters
        b_center = (self.b_max + self.b_min) / 2.0
        # adjust min and max based on integer value of naxis2
        self.b_max = b_center + (across_naxis / 2.0) * across_cdelt
        self.b_min = b_center - (across_naxis / 2.0) * across_cdelt

        across_coord = np.zeros(across_naxis)
        start = self.b_min + across_cdelt / 2.0
        for i in range(across_naxis):
            across_coord[i] = start
            start = start + across_cdelt

        if self.instrument == 'MIRI':
            self.crval1 = self.a_min
            self.xcoord = acoord
            self.ycoord = across_coord
            self.crval2 = self.b_min

        if self.instrument == 'NIRSPEC':
            self.crval2 = self.a_min
            self.ycoord = acoord
            self.xcoord = across_coord
            self.crval1 = self.b_min
# _______________________________________________________________________
# common to both MIRI and NIRSPEC
        self.crpix1 = 0.5
        self.crpix2 = 0.5

        # set up the lambda (z) coordinate of the cube
        range_lambda = self.lambda_max - self.lambda_min
        self.naxis3 = int(math.ceil(range_lambda / self.cdelt3))

        # adjust max based on integer value of naxis3
        lambda_center = (self.lambda_max + self.lambda_min) / 2.0

        self.lambda_min = lambda_center - (self.naxis3 / 2.0) * self.cdelt3
        self.lambda_max = lambda_center + (self.naxis3 / 2.0) * self.cdelt3

        self.lambda_max = self.lambda_min + (self.naxis3) * self.cdelt3
        self.zcoord = np.zeros(self.naxis3)
        zstart = self.lambda_min + self.cdelt3 / 2.0
        self.crval3 = zstart
        self.crpix3 = 1.0
        for i in range(self.naxis3):
            self.zcoord[i] = zstart
            zstart = zstart + self.cdelt3
# _______________________________________________________________________

    def print_cube_geometry(self):
        """Print out the general properties of the size of the IFU Cube
        """

        log.info('Cube Geometry:')
        if self.coord_system == 'internal_cal':
            log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(arcsec)  Min & Max (along slice, across slice)')
        else:
            log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(arcsec)  Min & Max (xi, eta arcsec)')
        log.info('Axis 1 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                 self.naxis1, self.crpix1, self.crval1, self.cdelt1,
                 self.a_min, self.a_max)
        log.info('Axis 2 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                 self.naxis2, self.crpix2, self.crval2, self.cdelt2,
                 self.b_min, self.b_max)
        if self.linear_wavelength:
            log.info('axis#  Naxis  CRPIX    CRVAL      CDELT(microns)  Min & Max (microns)')
            log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f %12.8f',
                     self.naxis3, self.crpix3, self.crval3, self.cdelt3,
                     self.lambda_min, self.lambda_max)

        if not self.linear_wavelength:
            log.info('Non-linear wavelength dimension; CDELT3 variable')
            log.info('axis#  Naxis  CRPIX    CRVAL     Min & Max (microns)')
            log.info('Axis 3 %5d  %5.2f %12.8f %12.8f %12.8f',
                     self.naxis3, self.crpix3, self.crval3,
                     self.wavelength_table[0], self.wavelength_table[self.naxis3 - 1])

        if self.rot_angle is not None:
            log.info('Rotation angle between Ra-Dec and Slicer-Plane %12.8f', self.rot_angle)

        if self.instrument == 'MIRI':
            # length of channel and subchannel are the same
            number_bands = len(self.list_par1)
            for i in range(number_bands):
                this_channel = self.list_par1[i]
                this_subchannel = self.list_par2[i]
                log.info(f'Cube covers channel, subchannel: {this_channel}, {this_subchannel}')
        elif self.instrument == 'NIRSPEC':
            # number of filters and gratings are the same
            number_bands = len(self.list_par1)
            for i in range(number_bands):
                this_fwa = self.list_par2[i]
                this_gwa = self.list_par1[i]
                log.info(f'Cube covers grating, filter: {this_gwa}, {this_fwa}')
# ________________________________________________________________________________

    def build_ifucube(self):
        """ Create the IFU cube

        1. Loop over every band contained in the IFU cube and read in the data
        associated with the band
        2. map_detector_to_output_frame: Maps the detector data to the cube output coordinate system
        3. For each mapped detector pixel the ifu cube spaxel located in the region of
        interest. There are two different routines to do this step, both of which use a c extension
        to combine the detector fluxes that fall within a region of influence from the spaxel center
        a. src/cube_match_sky: This routine uses the modified
        Shepard method to determining the weighting function, which weights the detector
        fluxes based on the distance between the detector center and spaxel center.
        b. src/cube_match_internal is only for single exposure, single band cubes and
        the ifucube in created in the detector plane. The weighting function is based on
        the overlap of between the detector pixel and spaxel. This method is simplified
        to determine the overlap in the along slice-wavelength plane.
        4. find_spaxel_flux: find the final flux associated with each spaxel
        5. setup_final_ifucube_model
        6. output_ifucube

        Returns
        -------
        Returns an ifu cube

        """

        self.output_name = self.define_cubename()
        total_num = self.naxis1 * self.naxis2 * self.naxis3

        self.spaxel_flux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_weight = np.zeros(total_num, dtype=np.float64)
        self.spaxel_var = np.zeros(total_num, dtype=np.float64)
        self.spaxel_iflux = np.zeros(total_num, dtype=np.float64)
        self.spaxel_dq = np.zeros(total_num, dtype=np.uint32)
        # ______________________________________________________________________________

        nxyplane = self.naxis1 * self.naxis2

        if self.spaxel_z == -1 and self.spaxel_x == -1 and self.spaxel_y == -1:
            debug_cube_index = -1

        elif(self.spaxel_z < 0 or self.spaxel_x < 0 or self.spaxel_y < 0):
            print('Incorrect input for Debug Spaxel values. Counting starts at 0')
            debug_cube_index = -1
            print(self.spaxel_z, self.spaxel_x, self.spaxel_y)
        else:
            spaxel_z = self.spaxel_z
            spaxel_x = self.spaxel_x
            spaxel_y = self.spaxel_y
            debug_cube_index = spaxel_z * (nxyplane) + spaxel_y * self.naxis1 + spaxel_x
            log.info(f"Printing debug information for cube spaxel:  {spaxel_x} {spaxel_y} {spaxel_z}")

        # ______________________________________________________________________________
        subtract_background = True

        # loop over every file that covers this channel/subchannel (MIRI) or
        # Grating/filter(NIRSPEC)
        # and map the detector pixels to the cube spaxel

        number_bands = len(self.list_par1)
        k = 0
        for ib in range(number_bands):
            this_par1 = self.list_par1[ib]
            this_par2 = self.list_par2[ib]
            for input in self.master_table.FileMap[self.instrument][this_par1][this_par2]:
                # ________________________________________________________________________________
                # loop over the files that cover the spectral range the cube is for

                input_model = datamodels.open(input)
                self.input_models_this_cube.append(input_model.copy())
                # set up input_model to be first file used to copy in basic header info
                # to ifucube meta data
                if ib == 0 and k == 0:
                    input_model_ref = input_model

                log.debug(f"Working on Band defined by: {this_par1} {this_par2}")
                # --------------------------------------------------------------------------------
                # POINTCLOUD used for skyalign and IFUalign
                # --------------------------------------------------------------------------------
                if self.interpolation in ['pointcloud', 'drizzle']:
                    pixelresult = self.map_detector_to_outputframe(this_par1,
                                                                   subtract_background,
                                                                   input_model)

                    coord1, coord2, corner_coord, wave, dwave, flux, err, slice_no, rois_pixel, \
                        roiw_pixel, weight_pixel, softrad_pixel, scalerad_pixel, \
                        x_det, y_det = pixelresult

                    # by default flag the dq plane based on the FOV of the detector projected to sky
                    flag_dq_plane = 1
                    if self.skip_dqflagging:
                        flag_dq_plane = 0

                    # check that there is valid data returned
                    # If all the data is flagged as DO_NOT_USE - not common- then log warning
                    build_cube = True
                    if wave is None:
                        log.warning(f'No valid data found on file {input_model.meta.filename}')
                        flag_dq_plane = 0
                        build_cube = False
                    # ______________________________________________________________________
                    # C extension setup
                    # ______________________________________________________________________
                    start_region = 0
                    end_region = 0

                    if self.instrument == 'MIRI':
                        instrument = 0
                        start_region = self.instrument_info.GetStartSlice(this_par1)
                        end_region = self.instrument_info.GetEndSlice(this_par1)

                    else:  # NIRSPEC
                        instrument = 1

                    result = None
                    weight_type = 0  # default to emsm instead of msm
                    if self.weighting == 'msm':
                        weight_type = 1

                    if self.interpolation == 'pointcloud' and build_cube:
                        roiw_ave = np.mean(roiw_pixel)
                        result = cube_wrapper(instrument, flag_dq_plane, weight_type, start_region, end_region,
                                              self.overlap_partial, self.overlap_full,
                                              self.xcoord, self.ycoord, self.zcoord,
                                              coord1, coord2, wave, flux, err, slice_no,
                                              rois_pixel, roiw_pixel, scalerad_pixel,
                                              weight_pixel, softrad_pixel,
                                              self.cdelt3_normal,
                                              roiw_ave, self.cdelt1, self.cdelt2)

                        spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq = result
                        self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                        self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                        self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                        self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                        spaxel_dq.astype(np.uint)
                        self.spaxel_dq = np.bitwise_or(self.spaxel_dq, spaxel_dq)
                        result = None
                        del result
                        del spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq
                    if self.weighting == 'drizzle' and build_cube:
                        cdelt3_mean = np.nanmean(self.cdelt3_normal)
                        xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner_coord
                        linear = 0
                        if self.linear_wavelength:
                            linear = 1
                        result = cube_wrapper_driz(instrument, flag_dq_plane,
                                                   start_region, end_region,
                                                   self.overlap_partial, self.overlap_full,
                                                   self.xcoord, self.ycoord, self.zcoord,
                                                   coord1, coord2, wave, flux, err, slice_no,
                                                   xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4,
                                                   dwave,
                                                   self.cdelt3_normal,
                                                   self.cdelt1, self.cdelt2, cdelt3_mean, linear,
                                                   x_det, y_det, debug_cube_index)

                        spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq = result
                        self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                        self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                        self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                        self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                        spaxel_dq.astype(np.uint)
                        self.spaxel_dq = np.bitwise_or(self.spaxel_dq, spaxel_dq)
                        result = None
                        del result
                        del spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, spaxel_dq
                # --------------------------------------------------------------------------------
                #                     # AREA - 2d method only works for single files local slicer plane (internal_cal)
                # --------------------------------------------------------------------------------
                elif self.interpolation == 'area':
                    # --------------------------------------------------------------------------------
                    # MIRI
                    # --------------------------------------------------------------------------------
                    if self.instrument == 'MIRI':
                        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                              'alpha_beta')
                        start_region = self.instrument_info.GetStartSlice(this_par1)
                        end_region = self.instrument_info.GetEndSlice(this_par1)
                        regions = list(range(start_region, end_region + 1))

                        for i in regions:
                            log.info('Working on Slice # %d', i)
                            y, x = (det2ab_transform.label_mapper.mapper == i).nonzero()

                            # getting pixel corner - ytop = y + 1 (routine fails for y = 1024)
                            index = np.where(y < 1023)
                            y = y[index]
                            x = x[index]
                            slice = i - start_region
                            result = cube_internal_cal.match_det2cube(self.instrument,
                                                                      x, y, slice,
                                                                      input_model,
                                                                      det2ab_transform,
                                                                      self.xcoord, self.zcoord,
                                                                      self.crval1, self.crval3,
                                                                      self.cdelt1, self.cdelt3,
                                                                      self.naxis1, self.naxis2)
                            spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux = result
                            self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                            self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                            self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                            self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                            result = None
                            del spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, result
                    # --------------------------------------------------------------------------------
                    # NIRSPEC
                    # --------------------------------------------------------------------------------
                    if self.instrument == 'NIRSPEC':
                        nslices = 30

                        slicemap = [15, 14, 16, 13, 17, 12, 18, 11, 19, 10,
                                    20, 9, 21, 8, 22, 7, 23, 6, 24, 5, 25,
                                    4, 26, 3, 27, 2, 28, 1, 29, 0]

                        for i in range(nslices):

                            slice_wcs = nirspec.nrs_wcs_set_input(input_model, i)
                            x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
                            detector2slicer = slice_wcs.get_transform('detector', 'slicer')

                            result = cube_internal_cal.match_det2cube(self.instrument,
                                                                      x, y, slicemap[i],
                                                                      input_model,
                                                                      detector2slicer,
                                                                      self.ycoord, self.zcoord,
                                                                      self.crval2, self.crval3,
                                                                      self.cdelt2, self.cdelt3,
                                                                      self.naxis1, self.naxis2)
                            spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux = result
                            self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                            self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                            self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                            self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                            result = None
                            del spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, result
                k = k + 1
                input_model.close()
                del input_model
            # _______________________________________________________________________
            # done looping over files

        self.find_spaxel_flux()
        self.set_final_dq_flags()

        # shove Flux and iflux in the  final IFU cube
        result = self.setup_final_ifucube_model(input_model_ref)
        return result

    # ********************************************************************************
    def build_ifucube_single(self):
        """ Build a set of single mode IFU cubes used for outlier detection
        and background matching

        Loop over every band contained in the IFU cube and read in the data
        associated with the band. Map each band to the output cube  coordinate
        system

        """
        # loop over input models
        single_ifucube_container = ModelContainer()

        weight_type = 0  # default to emsm instead of msm
        if self.weighting == 'msm':
            weight_type = 1
        number_bands = len(self.list_par1)
        this_par1 = self.list_par1[0]  # single IFUcube only have a single channel
        j = 0
        for i in range(number_bands):
            this_par2 = self.list_par2[i]
            nfiles = len(self.master_table.FileMap[self.instrument][this_par1][this_par2])
            # ________________________________________________________________________________
            # loop over the files that cover the spectral range the cube is for
            for k in range(nfiles):
                input_model = self.master_table.FileMap[self.instrument][this_par1][this_par2][k]
                self.input_models_this_cube.append(input_model)
                log.debug("Working on next Single IFU Cube = %i" % (j + 1))

                # for each new data model create a new spaxel
                total_num = self.naxis1 * self.naxis2 * self.naxis3
                self.spaxel_flux = np.zeros(total_num, dtype=np.float64)
                self.spaxel_weight = np.zeros(total_num, dtype=np.float64)
                self.spaxel_iflux = np.zeros(total_num)
                self.spaxel_dq = np.zeros(total_num, dtype=np.uint32)
                self.spaxel_var = np.zeros(total_num, dtype=np.float64)

                subtract_background = False

                pixelresult = self.map_detector_to_outputframe(this_par1,
                                                               subtract_background,
                                                               input_model)

                coord1, coord2, corner_coord, wave, dwave, flux, err, slice_no, \
                    rois_pixel, roiw_pixel, weight_pixel, \
                    softrad_pixel, scalerad_pixel = pixelresult

                build_cube = True
                if wave is None:  # there is no valid data on the detector. Pixels are flagged as DO_NOT_USE.
                    build_cube = False
                # the following values are not needed in cube_wrapper because the DQ plane is not being
                # filled in
                flag_dq_plane = 0
                start_region = 0
                end_region = 0
                roiw_ave = 0

                if self.instrument == 'MIRI':
                    instrument = 0
                else:  # NIRSPEC
                    instrument = 1

                result = None

                if self.interpolation == 'pointcloud' and build_cube:
                    roiw_ave = np.mean(roiw_pixel)
                    result = cube_wrapper(instrument, flag_dq_plane, weight_type, start_region, end_region,
                                          self.overlap_partial, self.overlap_full,
                                          self.xcoord, self.ycoord, self.zcoord,
                                          coord1, coord2, wave, flux, err, slice_no,
                                          rois_pixel, roiw_pixel, scalerad_pixel,
                                          weight_pixel, softrad_pixel,
                                          self.cdelt3_normal,
                                          roiw_ave, self.cdelt1, self.cdelt2)
                    spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, _ = result

                    self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                    self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                    self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                    self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                    result = None
                    del result, spaxel_flux, spaxel_var, spaxel_iflux

                if self.weighting == 'drizzle' and build_cube:
                    cdelt3_mean = np.nanmean(self.cdelt3_normal)
                    xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4 = corner_coord
                    linear = 0
                    if self.linear_wavelength:
                        linear = 1
                    result = cube_wrapper_driz(instrument, flag_dq_plane,
                                               start_region, end_region,
                                               self.overlap_partial, self.overlap_full,
                                               self.xcoord, self.ycoord, self.zcoord,
                                               coord1, coord2, wave, flux, err, slice_no,
                                               xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4,
                                               dwave,
                                               self.cdelt3_normal,
                                               self.cdelt1, self.cdelt2, cdelt3_mean, linear)

                    spaxel_flux, spaxel_weight, spaxel_var, spaxel_iflux, _ = result
                    self.spaxel_flux = self.spaxel_flux + np.asarray(spaxel_flux, np.float64)
                    self.spaxel_weight = self.spaxel_weight + np.asarray(spaxel_weight, np.float64)
                    self.spaxel_var = self.spaxel_var + np.asarray(spaxel_var, np.float64)
                    self.spaxel_iflux = self.spaxel_iflux + np.asarray(spaxel_iflux, np.float64)
                    result = None
                    del result, spaxel_flux, spaxel_var, spaxel_iflux
                # ______________________________________________________________________
                # shove Flux and iflux in the  final ifucube
                self.find_spaxel_flux()

                # determine Cube Spaxel flux
                status = 0
                result = self.setup_final_ifucube_model(input_model)
                ifucube_model, status = result

                single_ifucube_container.append(ifucube_model)
                if status != 0:
                    log.debug("Possible problem with single ifu cube, no valid data in cube")
                j = j + 1
        return single_ifucube_container

    # **************************************************************************
    def determine_cube_parameters_internal(self):
        """Determine the spatial and spectral ifu size for coord_system = internal_cal

        """

        # ____________________________________________________________
        # internal_cal is for only 1 file and weighting= area
        # no msm or emsm  information is needed
        par1 = self.list_par1[0]
        par2 = self.list_par2[0]

        a_scale, b_scale, w_scale = self.instrument_info.GetScale(par1,
                                                                  par2)
        self.spatial_size = a_scale
        if self.scalexy != 0:
            self.spatial_size = self.scalexy

        min_wave = self.instrument_info.GetWaveMin(par1, par2)
        max_wave = self.instrument_info.GetWaveMax(par1, par2)
        if self.wavemin is None:
            self.wavemin = min_wave
        else:
            self.wavemin = np.float64(self.wavemin)
            self.wavemin_user = True

        if self.wavemax is None:
            self.wavemax = max_wave
        else:
            self.wavemax = np.float64(self.wavemax)
            self.wavemax_user = True

        if self.scalew != 0:
            self.spectral_size = self.scalew
            self.linear_wavelength = True
        else:
            self.spectral_size = w_scale
            self.linear_wavelength = True

# **************************************************************************

    def determine_cube_parameters(self):
        """Determine the spatial and wavelength roi size to use for
        selecting point cloud elements around the spaxel centers.

        If the IFU cube covers more than 1 band - then use the rules to
        define the Spatial and Wavelength roi size to use for the cube
        Current Rule: using the minimum

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
        scalerad = np.zeros(number_bands)
        minwave = np.zeros(number_bands)
        maxwave = np.zeros(number_bands)

        # ____________________________________________________________
        for i in range(number_bands):
            if self.instrument == 'MIRI':
                par1 = self.list_par1[i]
                par2 = self.list_par2[i]
            elif self.instrument == 'NIRSPEC':
                par1 = self.list_par1[i]
                par2 = self.list_par2[i]

            a_scale, b_scale, w_scale = self.instrument_info.GetScale(par1,
                                                                      par2)

            spaxelsize[i] = a_scale
            spectralsize[i] = w_scale
            minwave[i] = self.instrument_info.GetWaveMin(par1, par2)
            maxwave[i] = self.instrument_info.GetWaveMax(par1, par2)

            # pull out the values from the cube pars reference file
            roiw[i] = self.instrument_info.GetWaveRoi(par1, par2)
            rois[i] = self.instrument_info.GetSpatialRoi(par1, par2)
            power[i] = self.instrument_info.GetMSMPower(par1, par2)
            softrad[i] = self.instrument_info.GetSoftRad(par1, par2)
            scalerad[i] = self.instrument_info.GetScaleRad(par1, par2)

        # Check the spatial size. If it is the same for the array set up the parameters
        all_same = np.all(spaxelsize == spaxelsize[0])

        if all_same:
            self.spatial_size = spaxelsize[0]
            spatial_roi = rois[0]
        # if it is not the same then use the minimum value
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
            self.wavemin_user = True

        if self.wavemax is None:
            self.wavemax = max_wave
        else:
            self.wavemax = np.float64(self.wavemax)
            self.wavemax_user = True

        # now check spectral step - this will determine
        # if the wavelength dimension is linear or not
        all_same_spectral = np.all(spectralsize == spectralsize[0])

        # check if scalew has been set - if yes then linear scale

        if self.scalew != 0:
            self.spectral_size = self.scalew
            self.linear_wavelength = True
            wave_roi = np.amin(roiw)
            weight_power = np.amin(power)
            self.soft_rad = np.amin(softrad)
            self.scalerad = np.amin(scalerad)

        # if all bands have the same spectral size then linear_wavelength
        elif all_same_spectral:
            self.spectral_size = spectralsize[0]
            wave_roi = roiw[0]
            weight_power = power[0]
            self.linear_wavelength = True  # added this 10/01/19
            self.soft_rad = softrad[0]
            self.scalerad = scalerad[0]
        else:
            self.linear_wavelength = False
            if self.instrument == 'MIRI':

                table = self.instrument_info.Get_multichannel_table(self.weighting)
                (table_wavelength, table_sroi,
                 table_wroi, table_power,
                 table_softrad, table_scalerad) = table

            # getting NIRSPEC Table Values
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
                    (table_wavelength, table_sroi,
                     table_wroi, table_power,
                     table_softrad, table_scalerad) = table
            # based on Min and Max wavelength - pull out the tables values that fall in this range
            # find the closest table entries to the self.wavemin and self.wavemax limits
            imin = (np.abs(table_wavelength - self.wavemin)).argmin()
            imax = (np.abs(table_wavelength - self.wavemax)).argmin()

            if imin > 1 and table_wavelength[imin] > self.wavemin:
                imin = imin - 1
            if (imax < (len(table_wavelength) - 1) and
                    self.wavemax > table_wavelength[imax]):
                imax = imax + 1

            num_table = imax - imin + 1
            self.scalerad_table = np.zeros(num_table)
            self.scalerad_table[:] = table_scalerad[imin:imax + 1]

            self.softrad_table = np.zeros(num_table)
            self.softrad_table[:] = table_softrad[imin:imax + 1]

            self.roiw_table = np.zeros(num_table)
            self.roiw_table[:] = table_wroi[imin:imax + 1]

            self.rois_table = np.zeros(num_table)
            self.rois_table[:] = table_sroi[imin:imax + 1]
            if self.num_files < 4:
                self.rois_table = self.rois_table * 1.5

            self.weight_power_table = np.zeros(num_table)
            self.weight_power_table[:] = table_power[imin:imax + 1]

            self.wavelength_table = np.zeros(num_table)
            self.wavelength_table[:] = table_wavelength[imin:imax + 1]

        # check if using default values from the table  (not user set)

        if self.rois == 0.0:
            self.rois = spatial_roi
            # not set by user but determined from tables
            # default rois in tables is designed with a 4 dither pattern
            # increase rois if less than 4 file

            if self.output_type == 'single' or self.num_files < 4:
                # We don't need to increase it if using 'emsm' weighting
                if self.weighting.lower() != 'emsm':
                    self.rois = self.rois * 1.5
                    log.info('Increasing spatial region of interest '
                             f'default value set for 4 dithers {self.rois}')

                # set wave_roi and  weight_power to same values if they are in  list
        if self.roiw == 0:
            self.roiw = wave_roi
        if self.weight_power == 0:
            self.weight_power = weight_power
        if self.scalexy != 0:
            self.spatial_size = self.scalexy

        # check on valid values

        found_error = False
        if self.linear_wavelength:
            # check we have valid data for key values
            if self.interpolation == 'pointcloud':
                if np.isnan(self.rois):
                    log.error('Spatial roi is nan, possible reference file value error')
                    found_error = True
                if np.isnan(self.roiw):
                    log.error('Spectral roi is nan, possible reference file value error')
                    found_error = True

            if self.weighting == 'msm':
                if np.isnan(self.weight_power):
                    log.error('Weight power is nan, possible reference file value error')
                    found_error = True
                if np.isnan(self.soft_rad):
                    log.error('Soft rad is nan, possible reference file value error')
                    found_error = True
            if self.weighting == 'emsm':
                if np.isnan(self.scalerad):
                    log.error('Scalerad is nan, possible reference file value error')
                    found_error = True
        else:
            if np.isnan(self.wavelength_table).all():
                log.error('Wavelength table contains all nans, possible reference file value error')
                found_error = True
            if self.interpolation == 'pointcloud':
                if np.isnan(self.rois_table).all():
                    log.error('Spatial roi table contains all nans, possible reference file value error')
                    found_error = True
                if np.isnan(self.roiw_table).all():
                    log.error('Spectral roi table contains all nans, possible reference file value error')
                    found_error = True

            if self.weighting == 'msm':
                if np.isnan(self.softrad_table).all():
                    log.error('Soft rad table contains all nans, possible reference file value error')
                    found_error = True
                if np.isnan(self.weight_power_table).all():
                    log.error('Weight power table contains all nans, possible reference file value error')
                    found_error = True
            if self.weighting == 'emsm':
                if np.isnan(self.scalerad_table).all():
                    log.error('Scalerad table contains all nans, possible reference file value error')
                    found_error = True

        if found_error:
            raise IncorrectParameter("An essential parameter is = nan, refer to apply error message")

        # catch where self.weight_power = nan weighting = msm written to header
        # TODO update writing to header scalerad if weighting = emsm

        if self.weight_power is not None:
            if np.isnan(self.weight_power):
                self.weight_power = None

        log.debug(f'spatial size {self.spatial_size}')
        log.debug(f'spectral size {self.spectral_size}')
        log.debug(f'spatial roi {self.rois}')
        log.debug(f'wave min and max {self.wavemin} {self.wavemax}')
        log.debug(f'linear wavelength {self.linear_wavelength}')
        log.debug(f'roiw {self.roiw}')
        log.debug(f'output_type {self.output_type}')
        if self.weighting == 'msm':
            log.debug(f'weight_power {self.weight_power}')
            log.debug(f'softrad {self.soft_rad}')
        if self.weighting == 'emsm':
            log.debug(f'scalerad {self.scalerad}')

# ******************************************************************************

    def setup_ifucube_wcs(self):
        """Function to determine the min and max coordinates of the spectral
        cube

        Loop over every datamodel contained in the cube and find the WCS
        of the output cube that contains all the data.

        Returns
        -------
        Footprint of cube: min and max of coordinates of cube.

        Notes
        -----
        If the coordinate system is internal_cal then min and max
        coordinates of along slice, across slice  and lambda (microns)

        For MIRI the units along/across slice dimension are arc seconds
        For NIRSPEC the units along/across slice dimension are meters

        If the coordinate system is skyalign/ifualign then the min and max of
        ra(degrees), dec (degrees) and lambda (microns) is returned.

        """
# _____________________________________________________________________________
        self.cdelt1 = self.spatial_size
        self.cdelt2 = self.spatial_size
        if self.linear_wavelength:
            self.cdelt3 = self.spectral_size

        parameter1 = self.list_par1
        parameter2 = self.list_par2
# ________________________________________________________________________________
# Define the rotation angle

# If coord_system = ifualign then the angle is between the ra-dec and alpha beta
# coord system using the first input model. Use first file in first band to set up rotation angle
# Compute the rotation angle between local IFU system  and RA-DEC

        if self.coord_system == 'ifualign':
            this_a = parameter1[0]  # 0 is first band - this_a is channel
            this_b = parameter2[0]  # 0 is first band - this_b is sub-channel
            log.info(f'Defining rotation between ra-dec and IFU plane using {this_a}, {this_b}')
            # first file for this band
            input_model = self.master_table.FileMap[self.instrument][this_a][this_b][0]

            if self.instrument == 'MIRI':
                xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_a)
                ysize = input_model.data.shape[0]
                y, x = np.mgrid[:ysize, xstart:xend]
                detector2alpha_beta = input_model.meta.wcs.get_transform('detector',
                                                                         'alpha_beta')
                alpha, beta, lam = detector2alpha_beta(x, y)
                lam_med = np.nanmedian(lam)
                # pick two alpha, beta values to determine rotation angle
                # values in arc seconds
                alpha_beta2world = input_model.meta.wcs.get_transform('alpha_beta',
                                                                      input_model.meta.wcs.output_frame.name)

                temp_ra1, temp_dec1, lam_temp = alpha_beta2world(0, 0, lam_med)
                temp_ra2, temp_dec2, lam_temp = alpha_beta2world(2, 0, lam_med)

            elif self.instrument == 'NIRSPEC':
                slice_wcs = nirspec.nrs_wcs_set_input(input_model, 0)
                x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
                detector2slicer = slice_wcs.get_transform('detector', 'slicer')
                across, along, lam = detector2slicer(x, y)  # lam ~0 for this transform
                lam_med = np.nanmedian(lam)

                # pick two along slice, across slice  values to determine rotation angle
                # values in meters
                slicer2world = slice_wcs.get_transform('slicer', slice_wcs.output_frame.name)
                temp_ra1, temp_dec1, lam_temp = slicer2world(0, 0, lam_med)
                temp_ra2, temp_dec2, lam_temp = slicer2world(0, 0.005, lam_med)
            # ________________________________________________________________________________
            # temp_dec1 is in degrees
            dra, ddec = (temp_ra2 - temp_ra1) * np.cos(temp_dec1 * np.pi / 180.0), (temp_dec2 - temp_dec1)
            self.rot_angle = 90 + np.arctan2(dra, ddec) * 180. / np.pi
            log.info(f'Rotation angle between ifu and sky: {self.rot_angle}')

        # If coord_system = iskyalign and the user provided a position angle. Define the rotation angle
        # to be the user provided value.

        if self.coord_system == 'skyalign' and self.cube_pa is not None:
            self.rot_angle = self.cube_pa
            log.info(f'Setting rotation angle between ifu and sky: {self.rot_angle}')
# ________________________________________________________________________________
# now loop over data and find min and max ranges data covers

        corner_a = []
        corner_b = []
        lambda_min = []
        lambda_max = []

        self.num_bands = len(self.list_par1)
        log.debug('Number of bands in cube: %i', self.num_bands)

        for i in range(self.num_bands):
            this_a = parameter1[i]
            this_b = parameter2[i]
            log.debug(f'Working on data from {this_a}, {this_b}')
            n = len(self.master_table.FileMap[self.instrument][this_a][this_b])
            log.debug('number of files %d', n)
            for k in range(n):
                lmin = 0.0
                lmax = 0.0

                input_file = self.master_table.FileMap[self.instrument][this_a][this_b][k]
                input_model = datamodels.open(input_file)

                # Find the footprint of the image
                spectral_found = hasattr(input_model.meta.wcsinfo, 'spectral_region')
                spatial_found = hasattr(input_model.meta.wcsinfo, 's_region')
                world = False
                if self.coord_system == 'skyalign' and self.cube_pa is None:
                    world = True

                # Do not use the default spatial or spectral region found in the wcs if
                # 1. instrument is MIRI and
                # 2. Output type is not multi and (not default calspec2) and
                # 3. Channel is 1 or 3 - channel with smaller FOV on detector
                if self.instrument == 'MIRI' and self.output_type != 'multi':
                    ch1 = '1'
                    ch3 = '3'
                    if ch1 in self.list_par1 or ch3 in self.list_par1:
                        spatial_found = False
                        spectral_found = False

                # If Moving Target data, then do not use S_REGION values.
                # The S_REGION values contain the footprint
                # on the sky of the original WCS.
                target_type = input_model.meta.target.type
                if target_type == 'MOVING':
                    spatial_found = False
                if spectral_found and spatial_found and world:
                    [lmin, lmax] = input_model.meta.wcsinfo.spectral_region
                    spatial_box = input_model.meta.wcsinfo.s_region
                    s = spatial_box.split(' ')
                    ca1 = float(s[3])
                    cb1 = float(s[4])
                    ca2 = float(s[5])
                    cb2 = float(s[6])
                    ca3 = float(s[7])
                    cb3 = float(s[8])
                    ca4 = float(s[9])
                    cb4 = float(s[10])
                else:
                    log.info('Mapping all pixels to output to determine IFU foot print')

                    if self.instrument == 'NIRSPEC':
                        ch_corners = cube_build_wcs_util.find_corners_NIRSPEC(
                            input_model,
                            self.instrument_info,
                            self.coord_system)
                        ca1, cb1, ca2, cb2, ca3, cb3, ca4, cb4, lmin, lmax = ch_corners
                    if self.instrument == 'MIRI':
                        ch_corners = cube_build_wcs_util.find_corners_MIRI(
                            input_model,
                            this_a,
                            self.instrument_info,
                            self.coord_system)

                        ca1, cb1, ca2, cb2, ca3, cb3, ca4, cb4, lmin, lmax = ch_corners

                # now append this model spatial and spectral corner
                corner_a.append(ca1)
                corner_a.append(ca2)
                corner_a.append(ca3)
                corner_a.append(ca4)

                corner_b.append(cb1)
                corner_b.append(cb2)
                corner_b.append(cb3)
                corner_b.append(cb4)

                lambda_min.append(lmin)
                lambda_max.append(lmax)
                input_model.close()
        # ________________________________________________________________________________
        # done looping over files determine final size of cube
        corner_a = np.array(corner_a)
        corner_a = wrap_ra(corner_a)

        final_a_min = min(corner_a)
        final_a_max = max(corner_a)
        final_b_min = min(corner_b)
        final_b_max = max(corner_b)

        log.debug(f'final a and b:{final_a_min, final_b_min, final_a_max, final_b_max}')
        log.debug(f' min and max wavelengths:   {min(lambda_min), max(lambda_max)}')
        # ______________________________________________________________________
        # the wavelength limits of cube are determined from 1. User or 2. cubepars
        # reference file (in the priority)
        final_lambda_min = self.wavemin
        final_lambda_max = self.wavemax

        log.debug(f' final min and max used in IFUcube:   {final_lambda_min, final_lambda_max}')

        if self.instrument == 'MIRI' and self.coord_system == 'internal_cal':
            #  we have a 1 to 1 mapping in y across slice  dimension
            nslice = self.instrument_info.GetNSlice(parameter1[0])
            log.info(f'Across slice scale {self.cdelt2}')
            self.cdelt2 = (final_b_max - final_b_min) / nslice

            final_b_max = final_b_min + (nslice) * self.cdelt2
            log.info('Changed the across slice scale dimension so we have 1-1 mapping between b and slice #')
            log.info(f'New across slice scale {self.cdelt2}')

        # ______________________________________________________________________
        if self.instrument == 'NIRSPEC' and self.coord_system == 'internal_cal':
            #  we have a 1 to 1 mapping in x - across slice  dimension.
            nslice = 30
            log.info(f'Across slice scale {self.cdelt1}')
            self.cdelt1 = (final_b_max - final_b_min) / nslice
            final_b_max = final_b_min + (nslice) * self.cdelt1
            log.info('Changed the across slice scale dimension so we have 1-1 mapping between b and slice #')
            log.info(f'New across slice Scale {self.cdelt1}')
            self.cdelt2 = self.cdelt1 / 2.0

        # ________________________________________________________________________________
        # Test that we have data (NIRSPEC NRS2 only has IFU data for 3 configurations)
        test_a = final_a_max - final_a_min
        test_b = final_b_max - final_b_min
        tolerance1 = 0.00001
        if test_a < tolerance1 or test_b < tolerance1:
            log.info(f'No Valid IFU slice data found {test_a} {test_b}')
        # ________________________________________________________________________________
        # Based on Scaling and Min and Max values determine naxis1, naxis2, naxis3
        # set cube CRVALs, CRPIXs

        if self.coord_system == 'skyalign' or self.coord_system == 'ifualign':
            self.set_geometry(corner_a, corner_b, final_lambda_min, final_lambda_max)
        else:
            self.set_geometryAB(corner_a, corner_b, final_lambda_min, final_lambda_max)

        self.print_cube_geometry()

    # **************************************************************************
    def map_detector_to_outputframe(self, this_par1,
                                    subtract_background,
                                    input_model):

        """Loop over a file and map the detector pixels to the output cube

        The output frame is on the SKY (ra-dec)

        Return the coordinates of all the detector pixel in the output frame.
        In addition, an array of pixel fluxes and weighing parameters are
        determined. The pixel flux and weighing parameters are used later in
        the process to find the final flux of a cube spaxel based on the pixel
        fluxes and pixel weighing parameters that fall within the roi of
        spaxel center

        Parameters
        ----------
        this_par1 : str
           for MIRI this is the channel # for NIRSPEC this is the grating name
           only need for MIRI to distinguish which channel on the detector we have
        subtract_background : boolean
           if TRUE then subtract the background found in the mrs_imatch step only
           needed for MIRI data
        input: datamodel
           input data model

        Returns
        -------
        coord1 : numpy.ndarray
           coordinate for axis1 in output cube for mapped pixel
        coord2: numpy.ndarray
           coordinate for axis2 in output cube for mapped pixel
        wave: numpy.ndarray
           wavelength associated with coord1,coord2
        flux: numpy.ndarray
           flux associated with coord1, coord2
        err: numpy.ndarray
           err associated with coord1, coord2
        rois_det: float
           spatial roi size to use
        roiw_det: numpy.ndarray
           spectral roi size associated with coord1,coord2
        weight_det : numpy.ndarray
            weighting parameter association with coord1,coord2
        softrad_det : numpy.ndarray
            weighting parameter association with coord1,coord2
        """
        # initialize alpha_det and beta_det to None. These are filled in
        # if the instrument is MIRI and the weighting is miripsf

        wave_all = None
        slice_no_all = None  # Slice number
        dwave_all = None
        corner_coord_all = None

        # initializa values to be returned to None
        dwave = None
        corner_coord = None
        coord1 = None
        coord2 = None
        wave = None
        flux = None
        err = None
        slice_no = None
        rois_det = None
        roiw_det = None
        weight_det = None
        softrad_det = None
        scalerad_det = None

        if self.instrument == 'MIRI':
            sky_result = self.map_miri_pixel_to_sky(input_model, this_par1, subtract_background)
            (x, y, ra, dec, wave_all, slice_no_all, dwave_all, corner_coord_all) = sky_result

        elif self.instrument == 'NIRSPEC':
            sky_result = self.map_nirspec_pixel_to_sky(input_model)
            (x, y, ra, dec, wave_all, slice_no_all, dwave_all, corner_coord_all) = sky_result

        # ______________________________________________________________________________
        # The following is for both MIRI and NIRSPEC

        flux_all = input_model.data[y, x]
        err_all = input_model.err[y, x]
        dq_all = input_model.dq[y, x]
        valid2 = np.isfinite(flux_all)

        x_all = x
        y_all = y
        # Pre-select only data within a given wavelength range
        # This range is defined to include all pixels for which the chosen wavelength region
        # of interest would have them fall within one of the cube spectral planes
        # Note that the cube lambda refer to the spaxel midpoints, so we must account for both
        # the spaxel width and the ROI size

        if self.linear_wavelength:
            min_wave_tolerance = self.zcoord[0] - self.roiw
            max_wave_tolerance = self.zcoord[-1] + self.roiw
        else:
            min_wave_tolerance = self.zcoord[0] - np.max(self.roiw_table)
            max_wave_tolerance = self.zcoord[-1] + np.max(self.roiw_table)
        if self.interpolation == 'drizzle':
            dmax = np.nanmax(dwave_all)

            if self.linear_wavelength:
                min_wave_tolerance = self.zcoord[0] - (self.cdelt3 + dmax)
                max_wave_tolerance = self.zcoord[-1] + (self.cdelt3 + dmax)
            else:
                min_wave_tolerance = self.zcoord[0] - self.cdelt3_normal[0]
                max_wave_tolerance = self.zcoord[-1] + self.cdelt3_normal[-1]
        valid_min = np.where(wave_all >= min_wave_tolerance)
        not_mapped_low = wave_all.size - len(valid_min[0])
        valid_max = np.where(wave_all <= max_wave_tolerance)
        not_mapped_high = wave_all.size - len(valid_max[0])
        if not_mapped_low > 0:
            log.info('# of detector pixels not mapped to output plane: '
                     f'{not_mapped_low} with wavelength below {min_wave_tolerance}')

        if not_mapped_high > 0:
            log.info('# of detector pixels not mapped to output plane: '
                     f'{not_mapped_high} with wavelength above {max_wave_tolerance}')

        # ______________________________________________________________________________
        # using the DQFlags from the input_image find pixels that should be excluded
        # from the cube mapping

        valid3 = np.logical_and(wave_all >= min_wave_tolerance,
                                wave_all <= max_wave_tolerance)

        # find the location of good data

        bad1 = np.bitwise_and(dq_all, dqflags.pixel['DO_NOT_USE']).astype(bool)
        bad2 = np.bitwise_and(dq_all, dqflags.pixel['NON_SCIENCE']).astype(bool)
        good_data = np.where(~bad1 & ~bad2 & valid2 & valid3)

        num_good = len(good_data[0])
        if num_good == 0:  # This can occcur if all the pixels on the detector are marked DO_NOT_USE.

            return coord1, coord2, corner_coord, wave, dwave, flux, err, \
                slice_no, rois_det, roiw_det, weight_det, \
                softrad_det, scalerad_det

        # good data holds the location of pixels we want to map to cube
        # define variables as numpy arrays (numba needs this defined)
        flux_all_good = flux_all[good_data]
        good_shape = flux_all_good.shape
        flux = np.zeros(good_shape, dtype=np.float64)
        err = np.zeros(good_shape, dtype=np.float64)
        coord1 = np.zeros(good_shape, dtype=np.float64)
        coord2 = np.zeros(good_shape, dtype=np.float64)
        wave = np.zeros(good_shape, dtype=np.float64)
        slice_no = np.zeros(good_shape)

        flux[:] = flux_all_good
        err[:] = err_all[good_data]
        wave[:] = wave_all[good_data]
        slice_no[:] = slice_no_all[good_data]
        x_all = x_all[good_data]
        y_all = y_all[good_data]

        log.debug(f'After removing pixels based on criteria min and max wave: {np.min(wave)}, {np.max(wave)}')

        # based on the wavelength define the sroi, wroi, weight_power and
        # softrad to use in matching detector to spaxel values
        rois_det = np.zeros(wave.shape)
        roiw_det = np.zeros(wave.shape)
        weight_det = np.zeros(wave.shape)
        softrad_det = np.zeros(wave.shape)
        scalerad_det = np.zeros(wave.shape)

        # ________________________________________________________________________
        # if working with msm or emsm need roi sizes and other parameters defined:
        if self.weighting == 'msm' or self.weighting == 'emsm':
            if self.linear_wavelength:
                rois_det[:] = self.rois
                roiw_det[:] = self.roiw
                weight_det[:] = self.weight_power
                softrad_det[:] = self.soft_rad
                scalerad_det[:] = self.scalerad
            else:
                # for each wavelength find the closest point in the self.wavelength_table

                for iw, w in enumerate(wave):
                    self.find_closest_wave(iw, w,
                                           self.wavelength_table,
                                           self.rois_table,
                                           self.roiw_table,
                                           self.softrad_table,
                                           self.weight_power_table,
                                           self.scalerad_table,
                                           rois_det,
                                           roiw_det,
                                           softrad_det,
                                           weight_det,
                                           scalerad_det)
        # ________________________________________________________________________
        ra_use = ra[good_data]
        dec_use = dec[good_data]
        coord1, coord2 = coord.radec2std(self.crval1,
                                         self.crval2,
                                         ra_use, dec_use,
                                         self.rot_angle)

        if self.interpolation == 'drizzle':
            dwave = np.zeros(good_shape)
            dwave[:] = dwave_all[good_data]
            ra1 = corner_coord_all[0]
            dec1 = corner_coord_all[1]
            ra2 = corner_coord_all[2]
            dec2 = corner_coord_all[3]
            ra3 = corner_coord_all[4]
            dec3 = corner_coord_all[5]
            ra4 = corner_coord_all[6]
            dec4 = corner_coord_all[7]
            ra1 = ra1[good_data]
            dec1 = dec1[good_data]
            ra2 = ra2[good_data]
            dec2 = dec2[good_data]
            ra3 = ra3[good_data]
            dec3 = dec3[good_data]
            ra4 = ra4[good_data]
            dec4 = dec4[good_data]

            xi1, eta1 = coord.radec2std(self.crval1,
                                        self.crval2,
                                        ra1, dec1,
                                        self.rot_angle)

            xi2, eta2 = coord.radec2std(self.crval1,
                                        self.crval2,
                                        ra2, dec2,
                                        self.rot_angle)
            xi3, eta3 = coord.radec2std(self.crval1,
                                        self.crval2,
                                        ra3, dec3,
                                        self.rot_angle)
            xi4, eta4 = coord.radec2std(self.crval1,
                                        self.crval2,
                                        ra4, dec4,
                                        self.rot_angle)

            corner_coord = [xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4]
        return coord1, coord2, corner_coord, wave, dwave, flux, err, \
            slice_no, rois_det, roiw_det, weight_det, \
            softrad_det, scalerad_det, x_all, y_all
    # ______________________________________________________________________

    def map_miri_pixel_to_sky(self, input_model, this_par1, subtract_background):
        """Loop over a file and map the detector pixels to the output cube
        The output frame is on the SKY (ra-dec)

        Return the coordinates of all the detector pixel in the output frame.

        Parameters
        ----------
        this_par1 : str
           for MIRI this is the channel # for NIRSPEC this is the grating name
           only need for MIRI to distinguish which channel on the detector we have
        subtract_background : boolean
           if TRUE then subtract the background found in the mrs_imatch step only
           needed for MIRI data
        input: datamodel
           input data model

        Returns
        -------
        x, y, ra, dec, lambda, slice_no  of valid slice pixels

        """
        wave = None
        slice_no = None  # Slice number
        dwave = None
        corner_coord = None

        # check if background sky matching as been done in mrs_imatch step
        # If it has not been subtracted and the background has not been
        # subtracted - subtract it.
        num_ch_bgk = len(input_model.meta.background.polynomial_info)
        if num_ch_bgk > 0 and subtract_background and input_model.meta.background.subtracted is False:
            for ich_num in range(num_ch_bgk):
                poly = input_model.meta.background.polynomial_info[ich_num]
                poly_ch = poly.channel
                if poly_ch == this_par1:
                    apply_background_2d(input_model, poly_ch, subtract=True)

        # find the slice number of each pixel and fill in slice_det
        ysize, xsize = input_model.data.shape
        slice_det = np.zeros((ysize, xsize), dtype=int)
        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                              'alpha_beta')
        start_region = self.instrument_info.GetStartSlice(this_par1)
        end_region = self.instrument_info.GetEndSlice(this_par1)
        regions = list(range(start_region, end_region + 1))
        for i in regions:
            ys, xs = (det2ab_transform.label_mapper.mapper == i).nonzero()
            xind = _toindex(xs)
            yind = _toindex(ys)
            xind = np.ndarray.flatten(xind)
            yind = np.ndarray.flatten(yind)
            slice_det[yind, xind] = i

        # define the x,y detector values of channel to be mapped to desired coordinate system
        xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_par1)
        y, x = np.mgrid[:ysize, xstart:xend]
        y = np.reshape(y, y.size)
        x = np.reshape(x, x.size)

        # if self.coord_system == 'skyalign' or self.coord_system == 'ifualign':
        ra, dec, wave = input_model.meta.wcs(x, y)
        valid1 = ~np.isnan(ra)
        ra = ra[valid1]
        dec = dec[valid1]
        wave = wave[valid1]
        x = x[valid1]
        y = y[valid1]

        xind = _toindex(x)
        yind = _toindex(y)
        xind = np.ndarray.flatten(xind)
        yind = np.ndarray.flatten(yind)
        slice_no = slice_det[yind, xind]

        if self.interpolation == 'drizzle':
            # Delta wavelengths
            _, _, wave1 = input_model.meta.wcs(x, y - 0.4999)
            _, _, wave2 = input_model.meta.wcs(x, y + 0.4999)
            dwave = np.abs(wave1 - wave2)

            # Pixel corners
            pixfrac = 1.0
            alpha1, beta, _ = input_model.meta.wcs.transform('detector', 'alpha_beta', x - 0.4999 * pixfrac, y)
            alpha2, _, _ = input_model.meta.wcs.transform('detector', 'alpha_beta', x + 0.4999 * pixfrac, y)
            # Find slice width
            allbetaval = np.unique(beta)
            dbeta = np.abs(allbetaval[1] - allbetaval[0])
            ra1, dec1, _ = input_model.meta.wcs.transform('alpha_beta',
                                                          input_model.meta.wcs.output_frame, alpha1,
                                                          beta - dbeta * pixfrac / 2., wave)
            ra2, dec2, _ = input_model.meta.wcs.transform('alpha_beta',
                                                          input_model.meta.wcs.output_frame, alpha1,
                                                          beta + dbeta * pixfrac / 2., wave)
            ra3, dec3, _ = input_model.meta.wcs.transform('alpha_beta',
                                                          input_model.meta.wcs.output_frame, alpha2,
                                                          beta + dbeta * pixfrac / 2., wave)
            ra4, dec4, _ = input_model.meta.wcs.transform('alpha_beta',
                                                          input_model.meta.wcs.output_frame, alpha2,
                                                          beta - dbeta * pixfrac / 2., wave)

            corner_coord = [ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4]

        sky_result = (x, y, ra, dec, wave, slice_no, dwave, corner_coord)
        return sky_result

    # ______________________________________________________________________
    def map_nirspec_pixel_to_sky(self, input_model):

        """Loop over a file and map the detector pixels to the output cube

        The output frame is on the SKY (ra-dec)
        Return the coordinates of all the detector pixel in the output frame.

        Parameters
        ----------
        input: datamodel
        input data model

        Returns
        -------
        x, y, ra, dec, lambda, slice_no

        """
        # initialize the ra,dec, and wavelength arrays
        # we will loop over slice_nos and fill in values
        # the flag_det will be set when a slice_no pixel is filled in
        # at the end we will use this flag to pull out valid data

        ysize, xsize = input_model.data.shape
        ra_det = np.zeros((ysize, xsize))
        dec_det = np.zeros((ysize, xsize))
        lam_det = np.zeros((ysize, xsize))
        flag_det = np.zeros((ysize, xsize))
        slice_det = np.zeros((ysize, xsize), dtype=int)
        dwave_det = np.zeros((ysize, xsize))
        ra1_det = np.zeros((ysize, xsize))
        ra2_det = np.zeros((ysize, xsize))
        ra3_det = np.zeros((ysize, xsize))
        ra4_det = np.zeros((ysize, xsize))
        dec1_det = np.zeros((ysize, xsize))
        dec2_det = np.zeros((ysize, xsize))
        dec3_det = np.zeros((ysize, xsize))
        dec4_det = np.zeros((ysize, xsize))

        pixfrac = 1.0

        # determine the slice width using slice 1 and 3
        slice_wcs1 = nirspec.nrs_wcs_set_input(input_model, 0)
        detector2slicer = slice_wcs1.get_transform('detector', 'slicer')
        x, y = wcstools.grid_from_bounding_box(slice_wcs1.bounding_box)
        across1, along1, _ = detector2slicer(x, y - 0.4999 * pixfrac)
        across1 = across1[~np.isnan(across1)]
        slice_loc1 = np.unique(across1)

        slice_wcs3 = nirspec.nrs_wcs_set_input(input_model, 2)
        detector2slicer = slice_wcs3.get_transform('detector', 'slicer')
        x, y = wcstools.grid_from_bounding_box(slice_wcs3.bounding_box)
        across3, along3, _ = detector2slicer(x, y - 0.4999 * pixfrac)
        across3 = across3[~np.isnan(across3)]
        slice_loc3 = np.unique(across3)

        across_width = abs(slice_loc1 - slice_loc3)
        # for NIRSPEC each file has 30 slices
        # wcs information access separately for each slice
        nslices = 30
        log.info("Mapping each NIRSpec slice to sky for input file: %s", input_model.meta.filename)

        for ii in range(nslices):
            slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)
            x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)
            ra, dec, lam = slice_wcs(x, y)

            # the slices are curved on detector so a rectangular region returns NaNs
            valid = ~np.isnan(lam)
            x = x[valid]
            y = y[valid]
            ra = ra[valid]
            dec = dec[valid]
            lam = lam[valid]

            if self.interpolation == 'drizzle':
                # Delta wavelengths
                _, _, wave1 = slice_wcs(x - 0.4999, y)
                _, _, wave2 = slice_wcs(x + 0.4999, y)
                dwave = np.abs(wave1 - wave2)

                # Pixel corners
                pixfrac = 1.0
                detector2slicer = slice_wcs.get_transform('detector', 'slicer')
                slicer2world = slice_wcs.get_transform('slicer', slice_wcs.output_frame)
                across1, along1, lam1 = detector2slicer(x, y - 0.49 * pixfrac)
                across2, along2, lam2 = detector2slicer(x, y + 0.49 * pixfrac)

                # Ensure that our ordering wraps around the footprint instead of crossing
                # footprint on a diagonal
                ra1, dec1, _ = slicer2world(across1 - across_width * pixfrac / 2, along1, lam1)
                ra2, dec2, _ = slicer2world(across1 + across_width * pixfrac / 2, along1, lam1)

                ra3, dec3, _ = slicer2world(across2 + across_width * pixfrac / 2, along2, lam2)
                ra4, dec4, _ = slicer2world(across2 - across_width * pixfrac / 2, along2, lam2)

                # near the slice boundaries the corners can become Nan - do not use pixels with
                # Nan corners
                valid1 = np.logical_and(~np.isnan(ra1), ~np.isnan(ra2))
                valid2 = np.logical_and(~np.isnan(ra3), ~np.isnan(ra4))
                final = np.where(np.logical_and(valid1, valid2))

                x = x[final]
                y = y[final]
                ra = ra[final]
                dec = dec[final]
                lam = lam[final]
                ra1 = ra1[final]
                dec1 = dec1[final]
                ra2 = ra2[final]
                dec2 = dec2[final]
                ra3 = ra3[final]
                dec3 = dec3[final]
                ra4 = ra4[final]
                dec4 = dec4[final]
                dwave = dwave[final]

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
                slice_det[yind, xind] = ii + 1

                # fill in corner values
                dwave_det[yind, xind] = dwave
                ra1_det[yind, xind] = ra1
                ra2_det[yind, xind] = ra2
                ra3_det[yind, xind] = ra3
                ra4_det[yind, xind] = ra4

                dec1_det[yind, xind] = dec1
                dec2_det[yind, xind] = dec2
                dec3_det[yind, xind] = dec3
                dec4_det[yind, xind] = dec4

            else:   # not drizzling
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
                slice_det[yind, xind] = ii + 1

        # after looping over slices  - pull out valid values

        valid_data = np.where(flag_det == 1)
        y, x = valid_data

        ra = ra_det[valid_data]
        dec = dec_det[valid_data]
        wave = lam_det[valid_data]
        slice_no = slice_det[valid_data]
        dwave = dwave_det[valid_data]
        ra1 = ra1_det[valid_data]
        ra2 = ra2_det[valid_data]
        ra3 = ra3_det[valid_data]
        ra4 = ra4_det[valid_data]
        dec1 = dec1_det[valid_data]
        dec2 = dec2_det[valid_data]
        dec3 = dec3_det[valid_data]
        dec4 = dec4_det[valid_data]

        corner_coord = [ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4]
        sky_result = (x, y, ra, dec, wave, slice_no, dwave, corner_coord)
        return sky_result

    # ********************************************************************************
    def find_closest_wave(self, iw, w,
                          wavelength_table,
                          rois_table,
                          roiw_table,
                          softrad_table,
                          weight_power_table,
                          scalerad_table,
                          rois_det,
                          roiw_det,
                          softrad_det,
                          weight_det,
                          scalerad_det):

        """ Given a specific wavelength, find the closest value in the wavelength_table

        """
        ifound = (np.abs(wavelength_table - w)).argmin()
        rois_det[iw] = rois_table[ifound]
        roiw_det[iw] = roiw_table[ifound]
        softrad_det[iw] = softrad_table[ifound]
        weight_det[iw] = weight_power_table[ifound]
        scalerad_det[iw] = scalerad_table[ifound]

    # ********************************************************************************
    def find_spaxel_flux(self):

        """Depending on the interpolation method, find the flux for each spaxel value
        """
        # currently these are the same but in the future there could be a difference in
        # how the spaxel flux is determined according to self.interpolation.
        if self.interpolation == 'area':
            good = self.spaxel_iflux > 0
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]
            self.spaxel_var[good] = self.spaxel_var[good] / (self.spaxel_weight[good] * self.spaxel_weight[good])
        elif self.interpolation == 'pointcloud' or self.interpolation == 'drizzle':
            # Don't apply any normalization if no points contributed to a spaxel (i.e., don't divide by zero)
            good = self.spaxel_iflux > 0

            # Normalize the weighted sum of pixel fluxes by the sum of the weights
            self.spaxel_flux[good] = self.spaxel_flux[good] / self.spaxel_weight[good]
            # Normalize the variance by the square of the weights
            self.spaxel_var[good] = self.spaxel_var[good] / (self.spaxel_weight[good] * self.spaxel_weight[good])

    # ********************************************************************************
    def set_final_dq_flags(self):

        """ Set up the final dq flags, Good data(0) , NON_SCIENCE or DO_NOT_USE
        """

        # An initial set of dq flags was set in overlap_fov_with_spaxel or
        # overlap_slice_with_spaxel. The initial dq dlags are defined in ifu_cube
        # class:
        # self.overlap_partial = 4  # intermediate flag
        # self.overlap_full  = 2    # intermediate flag
        # self.overlap_hole = dqflags.pixel['DO_NOT_USE']
        # self.overlap_no_coverage = dqflags.pixel['NON_SCIENCE'] (also bitwise and with
        # dqflags.pixel['DO_NOT_USE'] )

        # compare the weight plane and spaxel_dq. The initial spaxel_dq flagging
        # has too small a FOV in NIRSpec line mapping case.

        # flatten to match the size of spaxel_weight
        self.spaxel_dq = np.ndarray.flatten(self.spaxel_dq)

        # the fov is an underestimate. Check the spaxel_weight plane
        # if weight map > 0 then set spaxel_dq to overlap_partial
        under_data = self.spaxel_weight > 0
        self.spaxel_dq[under_data] = self.overlap_partial

        # convert all remaining spaxel_dq of 0 to NON_SCIENCE + DO_NOT_USE
        # these pixel should have no overlap with the data
        non_science = self.spaxel_dq == 0
        self.spaxel_dq[non_science] = np.bitwise_or(self.overlap_no_coverage,
                                                    dqflags.pixel['DO_NOT_USE'])

        # refine where good data should be
        ind_full = np.where(np.bitwise_and(self.spaxel_dq, self.overlap_full))
        ind_partial = np.where(np.bitwise_and(self.spaxel_dq, self.overlap_partial))

        self.spaxel_dq[ind_full] = 0
        self.spaxel_dq[ind_partial] = 0

        location_holes = np.where((self.spaxel_dq == 0) & (self.spaxel_weight == 0))
        self.spaxel_dq[location_holes] = self.overlap_hole

        # one last check. Remove pixels flagged as hole but have 1 adjacent spaxel
        # that has no coverage (NON_SCIENCE).  If NON_SCIENCE flag is next to pixel
        # flagged as hole then set the Hole flag to NON_SCIENCE
        spaxel_dq_temp = self.spaxel_dq
        nxy = self.naxis1 * self.naxis2
        index = np.where(self.spaxel_dq == self.overlap_hole)
        for i in range(len(index[0])):
            iwave = int(index[0][i] / nxy)
            rem = index[0][i] - iwave * nxy
            yrem = int(rem / self.naxis1)
            xrem = rem - yrem * self.naxis1

            found = 0
            ij = 0
            # do not allow holes to occur at the edge of IFU cube
            if (yrem == 0 or yrem == (self.naxis2 - 1) or
                    xrem == 0 or xrem == (self.naxis1 - 1)):
                spaxel_dq_temp[index[0][i]] = np.bitwise_or(self.overlap_no_coverage,
                                                            dqflags.pixel['DO_NOT_USE'])
                found = 1
            # flag as NON_SCIENCE instead of hole if left, right, top, bottom pixel
            # is NON_SCIENCE
            xcheck = np.zeros(4, dtype=int)
            ycheck = np.zeros(4, dtype=int)
            # left
            xcheck[0] = xrem - 1
            ycheck[0] = yrem
            # right
            xcheck[1] = xrem + 1
            ycheck[1] = yrem
            # bottom
            xcheck[2] = xrem
            ycheck[2] = yrem - 1
            # top
            xcheck[3] = xrem
            ycheck[3] = yrem + 1

            while ((ij < 4) and (found == 0)):
                if (xcheck[ij] > 0 and xcheck[ij] < self.naxis1 and
                        ycheck[ij] > 0 and ycheck[ij] < self.naxis2):
                    index_check = iwave * nxy + ycheck[ij] * self.naxis1 + xcheck[ij]
                    # If the nearby spaxel_dq contains overlap_no_coverage
                    # then unmark dq flag as hole. A hole has to have nearby
                    # pixels all in FOV.
                    check = (np.bitwise_and(self.spaxel_dq[index_check],
                                            self.overlap_no_coverage) == self.overlap_no_coverage)
                    if check:
                        spaxel_dq_temp[index[0][i]] = np.bitwise_or(self.overlap_no_coverage,
                                                                    dqflags.pixel['DO_NOT_USE'])
                        found = 1
                ij = ij + 1

        self.spaxel_dq = spaxel_dq_temp
        location_holes = np.where(self.spaxel_dq == self.overlap_hole)
        ave_holes = len(location_holes[0]) / self.naxis3

        if ave_holes < 1:
            log.info('Average # of holes/wavelength plane is < 1')
        else:
            log.info('Average # of holes/wavelength plane: %i', ave_holes)
        log.info('Total # of holes for IFU cube is : %i', len(location_holes[0]))

    # ********************************************************************************
    def setup_final_ifucube_model(self, model_ref):
        """ Set up the final meta WCS info of IFUCube along with other fits keywords

        return IFUCube model

        """
        status = 0
        # loop over the wavelength planes to confirm each plane has some data
        # for initial or final planes that do not have any data - eliminated them
        # from the IFUcube
        # Rearrange values from 1d vectors into 3d cubes

        flux = self.spaxel_flux.reshape((self.naxis3,
                                         self.naxis2, self.naxis1))
        wmap = self.spaxel_iflux.reshape((self.naxis3,
                                          self.naxis2, self.naxis1))

        var = self.spaxel_var.reshape((self.naxis3,
                                       self.naxis2, self.naxis1))
        dq = self.spaxel_dq.reshape((self.naxis3,
                                     self.naxis2, self.naxis1))

        # For MIRI MRS, apply a quality cut to help fix spectral tearing at the ends of each band.
        # This is largely taken care of by the WCS regions file, but there will still be 1-2 possibly
        # problematic planes at the end of each band in multi-band cubes.
        # Do this by looking for how many good spaxels there are at each wavelength and finding outliers
        # from the trend.
        if self.instrument == 'MIRI':
            nz = flux.shape[0]
            # Create a vector of the number of good spaxels at each wavelength
            ngood = np.zeros(nz)
            for zz in range(0, nz):
                dqvec = dq[zz, :, :].ravel()
                good = np.where(dqvec == 0)
                ngood[zz] = len(good[0])
            # Find where this vector is non-zero, and compute 1% threshold of those good values
            good = np.where(ngood > 0)
            if len(good[0]) > 0:
                pctile = np.percentile(ngood[good], 3)
                # Figure out where the number of good values were less than 75% of threshold,
                # and zero out those arrays.
                lowcov = (np.where((ngood > 0) & (ngood < 0.75 * pctile)))[0]
                nlowcov = len(lowcov)
                log.info('Number of spectral tear planes adjusted: %i', nlowcov)
                for zz in range(0, nlowcov):
                    flux[lowcov[zz], :, :] = 0
                    wmap[lowcov[zz], :, :] = 0
                    var[lowcov[zz], :, :] = 0
                    dq[lowcov[zz], :, :] = dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['NON_SCIENCE']

        # Set np.nan values wherever the DO_NOT_USE flag is set
        dnu = np.where((dq & dqflags.pixel['DO_NOT_USE']) != 0)
        flux[dnu] = np.nan
        var[dnu] = np.nan

        var = np.sqrt(var)
        if self.linear_wavelength:
            ifucube_model = datamodels.IFUCubeModel(data=flux, dq=dq,
                                                    err=var,
                                                    weightmap=wmap)
        else:
            wave = np.asarray(self.wavelength_table, dtype=np.float32)
            num = len(wave)
            alldata = np.array(
                [(wave[None].T, )],
                dtype=[('wavelength', '<f4', (num, 1))]
            )

            ifucube_model = datamodels.IFUCubeModel(data=flux, dq=dq,
                                                    err=var,
                                                    weightmap=wmap,
                                                    wavetable=alldata)

        ifucube_model.update(model_ref)
        ifucube_model.meta.filename = self.output_name

        # Call model_blender if there are multiple inputs
        if len(self.input_models) > 1:
            saved_model_type = ifucube_model.meta.model_type
            self.blend_output_metadata(ifucube_model)
            # Reset to original
            ifucube_model.meta.model_type = saved_model_type
# ______________________________________________________________________
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
# ______________________________________________________________________
# single files are created for a single band,

        if self.output_type == 'single':
            with datamodels.open(model_ref) as input:
                # define the cubename for each single
                filename = input.meta.filename
                indx = filename.rfind('.fits')
                self.output_name_base = filename[:indx]
                self.output_file = None
                newname = self.define_cubename()
                ifucube_model.meta.filename = newname

                if self.instrument == 'MIRI':
                    outchannel = self.list_par1[0]
                    outband = self.list_par2[0]
                    ifucube_model.meta.instrument.channel = outchannel
                    ifucube_model.meta.instrument.band = outband.upper()
                else:
                    outgrating = self.list_par1[0]
                    ifucube_model.meta.instrument.grating = outgrating.upper()
# ______________________________________________________________________
        ifucube_model.meta.wcsinfo.crval1 = self.crval1
        ifucube_model.meta.wcsinfo.crval2 = self.crval2
        ifucube_model.meta.wcsinfo.crpix1 = self.crpix1
        ifucube_model.meta.wcsinfo.crpix2 = self.crpix2

        ifucube_model.meta.wcsinfo.cdelt1 = self.cdelt1 / 3600.0
        ifucube_model.meta.wcsinfo.cdelt2 = self.cdelt2 / 3600.0
        # Now that we've got a pixel scale, set photometric area keywords
        ifucube_model.meta.photometry.pixelarea_arcsecsq = (
            self.cdelt1 * self.cdelt2)
        ifucube_model.meta.photometry.pixelarea_steradians = (
            ifucube_model.meta.photometry.pixelarea_arcsecsq * 2.3504e-11)
        if self.linear_wavelength:
            ifucube_model.meta.wcsinfo.crval3 = self.crval3
            ifucube_model.meta.wcsinfo.cdelt3 = self.cdelt3
            ifucube_model.meta.wcsinfo.ctype3 = 'WAVE'
            ifucube_model.meta.wcsinfo.crpix3 = self.crpix3
            ifucube_model.meta.ifu.roi_spatial = float(self.rois)
            ifucube_model.meta.ifu.roi_wave = float(self.roiw)
        else:
            ifucube_model.meta.wcsinfo.ctype3 = 'WAVE-TAB'
            ifucube_model.meta.wcsinfo.ps3_0 = 'WCS-TABLE'
            ifucube_model.meta.wcsinfo.ps3_1 = 'wavelength'
            ifucube_model.meta.wcsinfo.crval3 = 1.0
            ifucube_model.meta.wcsinfo.crpix3 = 1.0
            ifucube_model.meta.wcsinfo.cdelt3 = None
            ifucube_model.meta.ifu.roi_wave = np.mean(self.roiw_table)
            ifucube_model.wavedim = '(1,{:d})'.format(num)

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

        if self.rot_angle is None:
            self.rot_angle = 0.
        ifucube_model.meta.wcsinfo.pc1_1 = -np.cos(self.rot_angle * np.pi / 180.)
        ifucube_model.meta.wcsinfo.pc1_2 = np.sin(self.rot_angle * np.pi / 180.)
        ifucube_model.meta.wcsinfo.pc2_1 = np.sin(self.rot_angle * np.pi / 180.)
        ifucube_model.meta.wcsinfo.pc2_2 = np.cos(self.rot_angle * np.pi / 180.)

        ifucube_model.meta.ifu.flux_extension = 'SCI'
        ifucube_model.meta.ifu.error_extension = 'ERR'
        ifucube_model.meta.ifu.error_type = 'ERR'
        ifucube_model.meta.ifu.dq_extension = 'DQ'
        ifucube_model.meta.ifu.weighting = str(self.weighting)
        # weight_power is needed for single cubes. Linear Wavelengths
        # if non-linear wavelengths then this will be None
        ifucube_model.meta.ifu.weight_power = self.weight_power

        with datamodels.open(model_ref) as input:
            ifucube_model.meta.bunit_data = input.meta.bunit_data
            ifucube_model.meta.bunit_err = input.meta.bunit_err

        if self.interpolation == 'drizzle':
            # stick in values of 0, otherwise it is NaN and
            # fits file can not be written because these
            # values are defined in ifucube.schema.yaml
            ifucube_model.meta.ifu.weight_power = 0
            ifucube_model.meta.ifu.roi_wave = 0
            ifucube_model.meta.ifu.roi_spatial = 0
            ifucube_model.meta.ifu.weighting = str(self.interpolation)

        if self.coord_system == 'internal_cal':
            # stick in values of 0, otherwise it is NaN and
            # fits file can not be written because these
            # values are defined in ifucube.schema.yaml
            ifucube_model.meta.ifu.weight_power = 0
            ifucube_model.meta.ifu.roi_wave = 0
            ifucube_model.meta.ifu.roi_spatial = 0

            # uncorrect cdelt for degree conversion
            ifucube_model.meta.wcsinfo.cdelt1 *= 3600.0
            ifucube_model.meta.wcsinfo.cdelt2 *= 3600.0

            # correct "RA" axis orientation
            ifucube_model.meta.wcsinfo.pc1_1 *= -1.0

            if self.instrument == 'MIRI':
                ifucube_model.meta.wcsinfo.cunit1 = 'arcsec'
                ifucube_model.meta.wcsinfo.cunit2 = 'arcsec'
                ifucube_model.meta.wcsinfo.ctype1 = 'MRSALPHA'
                ifucube_model.meta.wcsinfo.ctype2 = 'MRSBETA'

            if self.instrument == 'NIRSPEC':
                ifucube_model.meta.wcsinfo.cunit1 = 'meter'
                ifucube_model.meta.wcsinfo.cunit2 = 'meter'
                ifucube_model.meta.wcsinfo.ctype1 = 'NRSSLICERX'
                ifucube_model.meta.wcsinfo.ctype2 = 'NRSSLICERY'

        # set WCS information
        wcsobj = pointing.create_fitswcs(ifucube_model)
        ifucube_model.meta.wcs = wcsobj
        ifucube_model.meta.wcs.bounding_box = ((0, self.naxis1 - 1),
                                               (0, self.naxis2 - 1),
                                               (0, self.naxis3 - 1))

        ifucube_model.meta.cal_step.cube_build = 'COMPLETE'
        # problem with cube_build - contains only 0 data
        if status == 1:
            ifucube_model.meta.cal_step.cube_build = 'SKIPPED'

        result = (ifucube_model, status)
        return result

    # ********************************************************************************
    def blend_output_metadata(self, IFUCube):
        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = IFUCube.meta.filename
        blendmeta.blendmodels(IFUCube, inputs=self.input_models_this_cube,
                              output=output_file)


class IncorrectInput(Exception):
    """ Raises an exception if input parameter, Interpolation, is set to area
    when more than one file is used to build the cube.
    """
    pass


class IncorrectParameter(Exception):
    """ Raises an exception if cube building  parameter is nan
    """
    pass
