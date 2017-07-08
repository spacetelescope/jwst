
# Routines used for building cubes
from __future__ import absolute_import, print_function
import sys
import time
import numpy as np
import math
import json
import logging

from astropy.io import fits
from astropy.modeling import models
#from gwcs.utils import _domain_to_bounds
from ..associations import Association
from .. import datamodels
from ..assign_wcs import nirspec
from ..assign_wcs import pointing
from . import cube_build_io_util
from . import cube_build_wcs_util
from . import file_table
from . import instrument_defaults
from . import spaxel
from . import cube_overlap
from . import cube_cloud
from . import data_types

from gwcs import wcstools


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class CubeData(object):
# CubeData - holds all the importatn informtion for IFU Cube Building:
# wcs, data, reference data

    def __init__(self, cube_type,
                 input_models,
                 input_filenames,
                 output_name_base,
                 data_type,
                 par_filename,
                 resol_filename,
                 **pars):

        self.cube_type = cube_type
        self.input_models = input_models
        self.input_filenames = input_filenames
        self.output_name_base = output_name_base
        self.data_type = data_type
        self.par_filename = par_filename
        self.resol_filename = resol_filename
        

        self.single = pars.get('single')
        self.channel = pars.get('channel')
        self.subchannel = pars.get('subchannel')
        self.grating = pars.get('grating')
        self.filter = pars.get('filter')
        self.scale1 = pars.get('scale1')
        self.scale2 = pars.get('scale2')
        self.scalew = pars.get('scalew')
        self.rois = pars.get('rois')
        self.roiw = pars.get('roiw')
        self.output_file = pars.get('output_file')
        self.interpolation = pars.get('interpolation')
        self.coord_system = pars.get('coord_system')
        self.offset_list = pars.get('offset_list')
        self.wavemin = pars.get('wavemin')
        self.wavemax = pars.get('wavemax')
        self.weighting = pars.get('weighting')
        self.weight_power = pars.get('weight_power')
        self.xdebug = pars.get('xdebug')
        self.ydebug = pars.get('ydebug')
        self.zdebug = pars.get('zdebug')
        self.debug_pixel = pars.get('debug_pixel')
        self.spaxel_debug = pars.get('spaxel_debug')

        self.ra_offset = []  # units arc seconds
        self.dec_offset = [] # units arc seconds
        self.detector = None
        self.instrument = None
        self.num_bands = 0
        self.band_channel = []
        self.band_subchannel = []
        self.band_filter = []
        self.band_grating = []
        self.num_bands = 0
        self.output_name = ''
        self.number_files = 0

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
    def setup(self):

        """
        Short Summary
        -------------
        Set up the IFU cube
        Read in the input_models and fill in the dictionary master_table that stores
        the files for each channel/subchannel or grating/filter

        if the channel/subchannel or grating/filter is not set then determine which
        ones are found in the data

        Read in necessary reference data:
        * ra dec offset list
        * cube parameter reference file
        * if miripsf weighting paramter is set then read in resolution file

        Parameters
        ----------
        instrument_info holds the defaults roi sizes  for each channel/subchannel (MIRI)
        or grating (NIRSPEC)

        Returns
        -------
        self with necessary files filled in
        """
#________________________________________________________________________________
# Check if there is an offset list (this ra,dec dither offset list will probably
# only be used in testing)

        if self.data_type == 'singleton':
            self.offset_list = 'NA'
        if self.offset_list != 'NA':
            log.info('Going to read in dither offset list')
            cube_build_io_util.read_offset_file(self)
#________________________________________________________________________________
# Read in the input data (association table or single file)
# Fill in MasterTable   based on Channel/Subchannel  or filter/grating
# Also if there is an Offset list - fill in MasterTable.FileOffset
#________________________________________________________________________________
        master_table = file_table.FileTable()
        instrument, detector = master_table.set_file_table(self.input_models,
                                                           self.input_filenames,
                                                           self.ra_offset,
                                                           self.dec_offset)
#________________________________________________________________________________
# find out how many files are in the association table or if it is an single file
# store the input_filenames and input_models
        num = 0
        num = len(self.input_filenames)
        self.number_files = num
        self.detector = detector
        self.instrument = instrument
#________________________________________________________________________________
    # Determine which channels/subchannels or filter/grating cubes will be
    # constructed from.
    # fills in band_channel, band_subchannel, band_grating, band_filer
#________________________________________________________________________________
        cube_build_io_util.determine_band_coverage(self, master_table)
#________________________________________________________________________________
    # check on interpolation = area and coord_system=alpha-beta types of cubes
    # if interpolation = area also checks that the use did not supply a scale2
    # values (beta dim)
#________________________________________________________________________________
        cube_build_io_util.check_cube_type(self)

        self.output_name = cube_build_io_util.update_output_name(self)
        if not self.single:
            log.info('Output Name %s',self.output_name)
#            log.info('Output Base %s ', self.output_name_base)
#________________________________________________________________________________
# InstrumentDefaults is an  dictionary that holds default parameters for
# difference instruments and for each band
#________________________________________________________________________________
        instrument_info = instrument_defaults.InstrumentInfo()
#--------------------------------------------------------------------------------
        # Load the parameter ref file data model
        # fill in the appropriate fields in InstrumentInfo
        # with the cube parameters
        log.info('Reading  cube parameter file %s', self.par_filename)
        cube_build_io_util.read_cubepars(self,instrument_info)
#--------------------------------------------------------------------------------
        # Load the miri resolution ref file data model -
        # fill in the appropriate fields in instrument_info
        # with the cube parameters
        if(self.weighting == 'miripsf'):
            log.info('Reading default MIRI cube resolution file %s', self.resol_filename)
            cube_build_io_util.read_resolution_file(self,instrument_info)
#________________________________________________________________________________
# get the ROI sizes
        self.instrument_info = instrument_info
        roi = CubeData.determine_roi_size(self)
        # if the user has not set the size of the ROI then use defaults in reference
        # parameter file

        if self.roiw == 0.0: self.roiw = roi[0]
        if self.rois == 0.0: self.rois = roi[1]
        if self.interpolation == 'pointcloud':
            log.info('Region of interest  %f %f',self.rois,self.roiw)

#________________________________________________________________________________
# Set up values to return and acess for other parts of cube_build

        self.master_table = master_table
        
        return self.output_file

#********************************************************************************

    def setup_wcs(self):

#********************************************************************************
        """
        Short Summary
        -------------
        Function to determine the min and max coordinates of the spectral
        cube,given channel & subchannel


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
        if self.cube_type == 'File' or self.cube_type == 'ASN' :
            log.info('Building Cube %s ', self.output_name)

        # Scale is 3 dimensions and is determined from values held in  instrument_info.GetScale
        scale = cube_build_wcs_util.determine_scale(self)
        self.Cdelt1 = scale[0]
        self.Cdelt2 = scale[1]
        self.Cdelt3 = scale[2]

        if self.instrument == 'MIRI':
            parameter1 = self.band_channel
            parameter2 = self.band_subchannel
        elif self.instrument == 'NIRSPEC':
            parameter1 = self.band_grating
            parameter2 = self.band_filter

        a_min = []
        a_max = []
        b_min = []
        b_max = []
        lambda_min = []
        lambda_max = []

        log.info('Number of bands in cube  %i',self.num_bands)

        for i in range(self.num_bands):
            this_a = parameter1[i]
            this_b = parameter2[i]
            log.debug('Working on data  from %s,%s',this_a,this_b)
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

                with datamodels.ImageModel(ifile) as input_model:
                    if self.instrument == 'NIRSPEC':
                        flag_data = 0
                        ch_footprint = cube_build_wcs_util.find_footprint_NIRSPEC(self,
                                                                              input_model,
                                                                              flag_data)
                        amin, amax, bmin, bmax, lmin, lmax = ch_footprint
#________________________________________________________________________________
                    if self.instrument == 'MIRI':
                        ch_footprint = cube_build_wcs_util.find_footprint_MIRI(self,
                                                                           input_model,
                                                                           this_a,
                                                                           self.instrument_info)
                        amin, amax, bmin, bmax, lmin, lmax = ch_footprint

# If a dither offset list exists then apply the dither offsets (offsets in arc seconds)

                    amin = amin - c1_offset/3600.0
                    amax = amax - c1_offset/3600.0

                    bmin = bmin - c2_offset/3600.0
                    bmax = bmax - c2_offset/3600.0

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

        if(self.wavemin != None and self.wavemin > final_lambda_min):
            final_lambda_min = self.wavemin
            log.info('Changed min wavelength of cube to %f ',final_lambda_min)

        if(self.wavemax != None and self.wavemax < final_lambda_max):
            final_lambda_max = self.wavemax
            log.info('Changed max wavelength of cube to %f ',final_lambda_max)
#________________________________________________________________________________
        if self.instrument =='MIRI' and self.coord_system=='alpha-beta':
        #  we have a 1 to 1 mapping in beta dimension.
            nslice = self.instrument_info.GetNSlice(parameter1[0])
            log.info('Beta Scale %f ',self.Cdelt2)
            self.Cdelt2 = (final_b_max - final_b_min)/nslice
        #print('remove********')
        #Cube.Cdelt2 = 0.17722104
            final_b_max = final_b_min + (nslice)*self.Cdelt2
        #print('remove********')
            log.info('Changed the Beta Scale dimension so we have 1 -1 mapping between beta and slice #')
            log.info('New Beta Scale %f ',self.Cdelt2)
#________________________________________________________________________________
# Test that we have data (NIRSPEC NRS2 only has IFU data for 3 configurations)
        test_a = final_a_max - final_a_min
        test_b = final_b_max - final_b_min
        test_w = final_lambda_max - final_lambda_min
        tolerance1 = 0.00001
        tolerance2 = 0.1
        if(test_a < tolerance1 or test_b < tolerance1 or test_w < tolerance2):
            log.info('No Valid IFU slice data found %f %f %f ',test_a,test_b,test_w)
#________________________________________________________________________________
        cube_footprint = (final_a_min, final_a_max, final_b_min, final_b_max,
                     final_lambda_min, final_lambda_max)
#________________________________________________________________________________
    # Based on Scaling and Min and Max values determine naxis1, naxis2, naxis3
    # set cube CRVALs, CRPIXs and xyz coords (center  x,y,z vector spaxel centers)

        if(self.coord_system == 'ra-dec'):
            cube_build_wcs_util.set_geometry(self,cube_footprint)
        else:
            cube_build_wcs_util.set_geometryAB(self,cube_footprint) # local coordinate system

        cube_build_wcs_util.print_cube_geometry(self)


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
        self.spaxel = CubeData.create_spaxel(self)
        # now need to loop over every file that covers this channel/subchannel (MIRI)
        # or Grating/filter(NIRSPEC)
        #and map the detector pixels to the cube spaxel.
        if(self.instrument == 'MIRI'):
            parameter1 = self.band_channel
            parameter2 = self.band_subchannel
        elif(self.instrument == 'NIRSPEC'):
            parameter1 = self.band_grating
            parameter2 = self.band_filter

        number_bands = len(parameter1)
        t0 = time.time()
        for i in range(number_bands):
            this_par1 = parameter1[i]
            this_par2 = parameter2[i]

            log.debug("Working on Band defined by:%s %s " ,this_par1,this_par2)
            CubeData.map_detector_to_spaxel(self,this_par1, this_par2,self.spaxel)

        t1 = time.time()
        log.info("Time Map All slices on Detector to Cube = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
# Mapped all data to cube or Point Cloud
# now determine Cube Spaxel flux

        t0 = time.time()
        CubeData.find_spaxel_flux(self, self.spaxel)

        t1 = time.time()
        log.info("Time find Cube Flux= %.1f.s" % (t1 - t0,))

        IFUCube = CubeData.setup_IFUCube(self,0)


#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube
        CubeData.update_IFUCube(self,IFUCube, self.spaxel)

        print('***** ',IFUCube.meta.wcsinfo.crval1,IFUCube.meta.wcsinfo.cravl2,IFUCube.meta.wcsinfo.crval3)
        print('***** ',IFUCube.meta.wcs.crval1,IFUCube.meta.wcs.cravl2,IFUCube.meta.wcs.crval3)
        return IFUCube

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


        single_IFUCube = datamodels.ModelContainer()
        n = len(self.input_models)
        this_par1 = self.band_channel[0] # only one channel is used in this approach
        this_par2 = None # not import for this type of mapping

        self.weighting =='standard'
        c1_offset = 0
        c2_offset = 0
        for j in range(n):
            t0 = time.time()
# for each new data model create a new spaxel
            spaxel = []
            spaxel = CubeData.create_spaxel(self)

            with datamodels.ImageModel(self.input_models[j]) as input_model:
#********************************************************************************
# pulled necessary routines from   CubeData.map_detector_to_spaxel
                if self.instrument == 'MIRI':
#________________________________________________________________________________
                    xstart, xend = self.instrument_info.GetMIRISliceEndPts(this_par1)
                    y, x = np.mgrid[:1024, xstart:xend]
                    y = np.reshape(y, y.size)
                    x = np.reshape(x, x.size)

                    cube_cloud.match_det2cube(self,input_model,
                                              x, y, j,
                                              this_par1,this_par2,
                                              spaxel,
                                              c1_offset, c2_offset)

                elif instrument == 'NIRSPEC':
                    # each file, detector has 30 slices - wcs information access seperately for each slice
                    start_slice = 0
                    end_slice = 29
                    nslices = end_slice - start_slice + 1
                    regions = list(range(start_slice, end_slice + 1))
                    for ii in regions:
                        t0a = time.time()
                        slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)
                        yrange = slice_wcs.bounding_box[1][0],slice_wcs.bounding_box[1][1]
                        xrange = slice_wcs.bounding_box[0][0],slice_wcs.bounding_box[0][1]
                        x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box)


                        cube_cloud.match_det2cube(self,input_model,
                                                  x, y, ii,
                                                  this_par1,this_par2,
                                                  spaxel,
                                                  c1_offset, c2_offset)

                        t1a = time.time()
                        log.debug("Time Match one NIRSPEC slice  to IFUCube = %.1f.s" % (t1a - t0a,))
#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube
            CubeData.find_spaxel_flux(self, spaxel)
# now determine Cube Spaxel flux
            IFUCube = CubeData.setup_IFUCube(self,j)

            CubeData.update_IFUCube(self,IFUCube, spaxel)
            print('***** ',IFUCube.meta.wcsinfo.crval1,IFUCube.meta.wcsinfo.crval2,
                  IFUCube.meta.wcsinfo.crval3)

            print('wcs ',IFUCube.meta.wcs(1,1,1))
            

            t1 = time.time()
            log.info("Time Create Single IFUcube  = %.1f.s" % (t1 - t0,))
#            print('build_ifucube_single:',IFUCube.meta.filename)
#_______________________________________________________________________
            single_IFUCube.append(IFUCube)
            del spaxel[:]
        return single_IFUCube

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


        total_num = self.naxis1*self.naxis2*self.naxis3

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
            number_bands = len(self.band_channel)
            min_s = 1000.00
            min_w = 1000.00

            for i in range(number_bands):
                this_channel = self.band_channel[i]
                this_sub = self.band_subchannel[i]
                wroi = self.instrument_info.GetWaveRoi(this_channel,this_sub)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_channel,this_sub)
                if sroi < min_s:
                    min_s = sroi
            roi = [min_w, min_s]

        elif self.instrument == 'NIRSPEC':
            number_gratings = len(self.band_grating)
            min_s = 1000.00
            min_w = 1000.00

            for i in range(number_gratings):
                this_gwa = self.band_grating[i]
                wroi = self.instrument_info.GetWaveRoi(this_gwa)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_gwa)
                if sroi < min_s:
                    min_s = sroi
            roi = [min_w, min_s]
        return roi

#********************************************************************************
    def map_detector_to_spaxel(self,this_par1, this_par2,spaxel):
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

        instrument  = self.instrument
        nfiles = len(self.master_table.FileMap[instrument][this_par1][this_par2])
        log.debug('Number of files in cube %i', nfiles)

    # loop over the files that cover the spectral range the cube is for

        for k in range(nfiles):
            ifile = self.master_table.FileMap[instrument][this_par1][this_par2][k]
            ioffset = len(self.master_table.FileOffset[this_par1][this_par2]['C1'])

            c1_offset = 0.0
            c2_offset = 0.0
        # c1_offset and c2_offset are the dither offset sets (in arc seconds)
        # by default these are zer0. The user has to supply these
            if ioffset == nfiles:
                c1_offset = self.master_table.FileOffset[this_par1][this_par2]['C1'][k]
                c2_offset = self.master_table.FileOffset[this_par1][this_par2]['C2'][k]
# Open the input data model
            with datamodels.ImageModel(ifile) as input_model:
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
                        cube_cloud.match_det2cube(self,input_model,
                                            x, y, k,
                                            this_par1,this_par2,
                                            spaxel,
                                            c1_offset, c2_offset)


                        t1 = time.time()
                        log.debug("Time Match one Channel from 1 file  to IFUCube = %.1f.s"
                                  % (t1 - t0,))
#________________________________________________________________________________
#2D area method - only works for single files and coord_system = 'alpha-beta'
                    if self.interpolation == 'area':
                        det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                              'alpha_beta')


                        start_region = self.instrument_info.GetStartSlice(this_par1)
                        end_region = self.instrument_info.GetEndSlice(this_par1)
                        regions = list(range(start_region, end_region + 1))

                    #xtest = 28.310396-1 # test pixel to compare with Distortion doc
                    #ytest = 512.0-1     # test pixel to compare with Distortion doc
                    #coord1_test,coord2_test,lam_test = det2ab_shift(xtest,ytest)
                    #print('test values',xtest+1,ytest+1,coord1_test,coord2_test,lam_test)

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
                    start_slice = 0
                    end_slice = 29
                    nslices = end_slice - start_slice + 1
                    regions = list(range(start_slice, end_slice + 1))
                    log.info("Mapping each NIRSPEC slice to sky, this takes a while for NIRSPEC data")
                    for i in regions:
#                    print('on region ',i)
                        slice_wcs = nirspec.nrs_wcs_set_input(input_model, i)
                        yrange = slice_wcs.bounding_box[1][0],slice_wcs.bounding_box[1][1]
                        xrange = slice_wcs.bounding_box[0][0],slice_wcs.bounding_box[0][1]


                        x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box,
                                                              step=(1,1), center=True)
                        t0 = time.time()
                        cube_cloud.match_det2cube(self,input_model,
                                                  x, y, i,
                                                  this_par1,this_par2,
                                                  spaxel,
                                                  c1_offset, c2_offset)


                        t1 = time.time()
                        log.debug("Time Match one NIRSPEC slice  to IFUCube = %.1f.s" % (t1 - t0,))
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
                    spaxel[i].flux = spaxel[i].flux/spaxel[i].flux_weight

        elif self.interpolation == 'pointcloud':
            icube = 0
            t0 = time.time()
            for iz, z in enumerate(self.zcoord):
                for iy, y in enumerate(self.ycoord):
                    for ix, x in enumerate(self.xcoord):

                        if(spaxel[icube].iflux > 0):
                            spaxel[icube].flux = spaxel[icube].flux/spaxel[icube].flux_weight

                            if(self.debug_pixel == 1 and self.xdebug == ix and
                               self.ydebug == iy and self.zdebug == iz ):

                                log.debug('For spaxel %d %d %d final flux %f '
                                          %(self.xdebug+1,self.ydebug+1,
                                            self.zdebug+1,spaxel[icube].flux))
                                self.spaxel_debug.write('For spaxel %d %d %d, final flux %f '
                                                        %(self.xdebug+1,self.ydebug+1,
                                                          self.zdebug+1,spaxel[icube].flux) +' \n')
                        icube = icube + 1
            t1 = time.time()
            log.info("Time to interpolate at spaxel values = %.1f.s" % (t1 - t0,))

#********************************************************************************
    def setup_IFUCube(self,j):

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

        IFUCube = datamodels.IFUCubeModel(data=data, dq=dq_cube, err=err_cube, weightmap=idata)
        IFUCube.update(self.input_models[j])

        IFUCube.meta.filename = self.output_name
        if self.single:
            with datamodels.open(self.input_models[j]) as input:
                # makingf fileanme = org gives a error later when past
                # back to model container - do we want to define
                # a new KEYWORD - filename_org ?
                #IFUCube.meta.filename = input.meta.filename

                filename = self.input_filenames[j]
                indx = filename.rfind('.fits')
                self.output_name_base = filename[:indx]
                self.output_file = None
                newname  = cube_build_io_util.update_output_name(self)
                IFUCube.meta.filename = newname
                IFUCube.meta.instrument.channel = self.band_channel[0] 

        IFUCube.meta.wcsinfo.crval1 = self.Crval1
        IFUCube.meta.wcsinfo.crval2 = self.Crval2
        IFUCube.meta.wcsinfo.crval3 = self.Crval3
        IFUCube.meta.wcsinfo.crpix1 = self.Crpix1
        IFUCube.meta.wcsinfo.crpix2 = self.Crpix2
        IFUCube.meta.wcsinfo.crpix3 = self.Crpix3
        IFUCube.meta.wcsinfo.cdelt1 = self.Cdelt1/3600.0
        IFUCube.meta.wcsinfo.cdelt2 = self.Cdelt2/3600.0
        IFUCube.meta.wcsinfo.cdelt3 = self.Cdelt3

        IFUCube.meta.wcsinfo.ctype1 = 'RA---TAN'
        IFUCube.meta.wcsinfo.ctype2 = 'DEC--TAN'
        IFUCube.meta.wcsinfo.cunit1 = 'deg'
        IFUCube.meta.wcsinfo.cunit2 = 'deg'

        IFUCube.meta.wcsinfo.ctype3 = 'WAVE'
        IFUCube.meta.wcsinfo.cunit3 = 'um'
        IFUCube.meta.wcsinfo.wcsaxes = 3

        IFUCube.meta.flux_extension = 'SCI'
        IFUCube.meta.error_extension = 'ERR'
        IFUCube.meta.dq_extension = 'DQ'
#        IFUCube.meta.data_model_type = 'IFUCubeModel'
        IFUCube.error_type = 'ERR'

        if(self.coord_system == 'alpha-beta'):
            IFUCube.meta.wcsinfo.ctype1 = 'ALPHA'
            IFUCube.meta.wcsinfo.ctype2 = 'BETA'
            IFUCube.meta.wcsinfo.cunit1 = 'arcsec'
            IFUCube.meta.wcsinfo.cunit2 = 'arcsec'

        print('***** ',IFUCube.meta.wcsinfo.crval1,IFUCube.meta.wcsinfo.crval2,
                  IFUCube.meta.wcsinfo.crval3)

        wcsobj = pointing.create_fitswcs(IFUCube)
        IFUCube.meta.wcs = wcsobj
        IFUCube.meta.wcs.bounding_box = ((0,naxis1-1),(0,naxis2-1),(0,naxis3-1))
        return IFUCube

#********************************************************************************

    def update_IFUCube(self,IFUCube, spaxel):

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


        temp_flux =np.reshape(np.array([s.flux for s in spaxel]),
                          [self.naxis3,self.naxis2,self.naxis1])
        temp_wmap =np.reshape(np.array([s.iflux for s in spaxel]),
                          [self.naxis3,self.naxis2,self.naxis1])


        IFUCube.data = temp_flux
        IFUCube.weightmap = temp_wmap

        IFUCube.meta.cal_step.cube_build = 'COMPLETE'
#    icube = 0
#    for z in range(Cube.naxis3):
#        for y in range(Cube.naxis2):
#            for x in range(Cube.naxis1):
#                IFUCube.data[z, y, x] = spaxel[icube].flux
#                IFUCube.weightmap[z, y, x] = len(spaxel[icube].ipointcloud)
#                icube = icube + 1


       # result = IFUCube.copy()
        #return result

#********************************************************************************
