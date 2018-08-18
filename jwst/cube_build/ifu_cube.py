
# Routines used for building cubes
import sys
import time
import numpy as np
import logging
from ..model_blender import blendmeta
from .. import datamodels
from ..assign_wcs import nirspec
from ..assign_wcs import pointing
from . import cube_build_wcs_util
from . import spaxel
from . import cube_overlap
from . import cube_cloud
from gwcs import wcstools

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class IFUCubeData():
# CubeData - holds all the importatn informtion for IFU Cube Building:
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

        self.instrument  = instrument
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
# first define the number of file names that will be used to construct this cube
# do some checks on the IFUCube to be made
# find the ROI size

    def setup_cube(self):
        num1= len(self.list_par1)
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
            if(self.scale2 !=0):
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
               log.info('Increasing spatial region of interest' + \
                            ' default value set for 4 dithers %f', self.rois)
        if self.interpolation == 'pointcloud':
            log.info('Region of interest spatial, wavelength  %f %f',self.rois,self.roiw)

#________________________________________________________________________________

# update the output name
    def define_cubename(self):

#        print ('ifu_cube:define_cubename basename ',self.output_name_base)

        if self.pipeline == 2:
            newname  = self.output_name_base + '_s3d.fits'
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
                        if i < number_channels-1: ch_name = ch_name + '-'

                subchannels = list(set(self.list_par2))
                number_subchannels = len(subchannels)
                b_name = ''
                for i in range(number_subchannels):
                    b_name = b_name + subchannels[i]
                    if i > 1 : b_name = b_name + '-'
                b_name  = b_name.lower()
                newname = self.output_name_base + ch_name+'-'+ b_name + '_s3d.fits'
                if self.coord_system == 'alpha-beta':
                    newname = self.output_name_base + ch_name+'-'+ b_name + '_ab_s3d.fits'
                if self.output_type == 'single':
                    newname = self.output_name_base + ch_name+'-'+ b_name + '_single_s3d.fits'
#________________________________________________________________________________
            elif self.instrument == 'NIRSPEC':
                fg_name = '_'

                for i in range( len(self.list_par1)):
                    fg_name = fg_name + self.list_par1[i] + '-'+ self.list_par2[i]
                    if(i < self.num_bands -1):
                        fg_name = fg_name + '-'
                fg_name = fg_name.lower()
                newname = self.output_name_base + fg_name+ '_s3d.fits'
                if self.output_type == 'single':
                    newname = self.output_name_base + fg_name+ 'single_s3d.fits'
#________________________________________________________________________________
        if self.output_type != 'single':
            log.info('Output Name %s',newname)

#        print('*** newname ****',newname)
        return newname
#********************************************************************************
    class IncorrectInput(Exception):
        pass
    class AreaInterpolation(Exception):
        pass
#********************************************************************************
    def setup_ifucube_wcs(self):

        cube_build_wcs_util.setup_wcs(self)

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
        self.find_output_type()

        self.spaxel = self.create_spaxel()

        # now need to loop over every file that covers this channel/subchannel (MIRI)
        # or Grating/filter(NIRSPEC)
        #and map the detector pixels to the cube spaxel.

        number_bands = len(self.list_par1)
        t0 = time.time()
        for i in range(number_bands):
            this_par1 = self.list_par1[i]
            this_par2 = self.list_par2[i]

            log.debug("Working on Band defined by:%s %s " ,this_par1,this_par2)
            self.map_detector_to_spaxel(this_par1, this_par2,self.spaxel)

        t1 = time.time()
        log.info("Time Map All slices on Detector to Cube = %.1f.s" % (t1 - t0,))
#_______________________________________________________________________
# Mapped all data to cube or Point Cloud
# now determine Cube Spaxel flux

        t0 = time.time()
        self.find_spaxel_flux(self.spaxel)

        t1 = time.time()
        log.info("Time to find Cube Flux= %.1f.s" % (t1 - t0,))

        IFUCube = self.setup_IFUCube(0)
#_______________________________________________________________________
# shove Flux and iflux in the  final IFU cube
        self.update_IFUCube(IFUCube, self.spaxel)
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
        log.info("Number of Single IFU cubes creating  = %i" % n)
        this_par1 = self.list_par1[0] # only one channel is used in this approach
        this_par2 = None # not important for this type of mapping

        self.weighting =='msm'
        c1_offset = 0
        c2_offset = 0
        for j in range(n):
            log.info("Working on next Single IFU Cube  = %i" %(j+1))
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


                    cube_cloud.match_det2cube(self,input_model,
                                              x, y, j,
                                              this_par1,this_par2,
                                              spaxel,
                                              c1_offset, c2_offset)

                elif self.instrument == 'NIRSPEC':
                    # each file, detector has 30 slices - wcs information access seperately for each slice

                    nslices = 30 
                    for ii in range(nslices):
                        t0a = time.time()
                        slice_wcs = nirspec.nrs_wcs_set_input(input_model, ii)
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

            self.find_spaxel_flux(spaxel)
# now determine Cube Spaxel flux

            IFUCube = self.setup_IFUCube(j)
            self.update_IFUCube(IFUCube, spaxel)
            t1 = time.time()
            log.info("Time Create Single IFUcube  = %.1f.s" % (t1 - t0,))
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
            number_bands = len(self.list_par1)
            min_s = 1000.00
            min_w = 1000.00

            for i in range(number_bands):
                this_channel = self.list_par1[i]
                this_sub = self.list_par2[i]
                wroi = self.instrument_info.GetWaveRoi(this_channel,this_sub)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_channel,this_sub)
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
#                print('Grating and Filter',this_gwa,this_filter)

                wroi = self.instrument_info.GetWaveRoi(this_gwa,this_filter)
                if wroi < min_w:
                    min_w = wroi
                sroi = self.instrument_info.GetSpatialRoi(this_gwa,this_filter)
                if sroi < min_s:
                    min_s = sroi
            roi = [min_w, min_s]
        return roi

#********************************************************************************
    def map_detector_to_spaxel(self,this_par1, this_par2,spaxel):
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

        instrument  = self.instrument
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
                do_background_subtraction = False
                num_ch_bgk = len(input_model.meta.background.polynomial_info)

                if(num_ch_bgk> 0):

                    do_background_subtraction = True
                    for ich_num in range(num_ch_bgk):
                        poly = input_model.meta.background.polynomial_info[ich_num]
                        poly_ch = poly.channel
                        if(poly_ch == this_par1):
                            apply_background_2d(input_model,poly_ch,subtract=True)

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
#                        slice_wcs = nirspec.nrs_wcs_set_input(input_model, i)
#                        x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box,
#                                                              step=(1,1), center=True)
                        

                        t0 = time.time()
                        x = 0
                        y = 0 
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
        
        # Call model_blender if there are multiple inputs
        if len(self.input_models) > 1:
            saved_model_type = IFUCube.meta.model_type
            self.blend_output_metadata(IFUCube)
            IFUCube.meta.model_type = saved_model_type  # Reset to original

#______________________________________________________________________
        if self.output_type == 'single':
            with datamodels.open(self.input_models[j]) as input:
                # define the cubename for each single
                filename = self.input_filenames[j]
                indx = filename.rfind('.fits')
                self.output_name_base = filename[:indx]
                self.output_file = None
                newname  = self.define_cubename()
                IFUCube.meta.filename = newname

#______________________________________________________________________
# fill in Channel, Band for MIRI
        if self.instrument == 'MIRI':
            # fill in Channel output meta data
            num_ch = len(self.list_par1)
            IFUCube.meta.instrument.channel = self.list_par1[0]
            num_ch = len(self.list_par1)
            for m in range (1, num_ch):
                IFUCube.meta.instrument.channel =  IFUCube.meta.instrument.channel + str(self.list_par1[m])

#______________________________________________________________________
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
        IFUCube.meta.wcsinfo.pc1_1 = -1
        IFUCube.meta.wcsinfo.pc1_2 = 0
        IFUCube.meta.wcsinfo.pc1_3 = 0

        IFUCube.meta.wcsinfo.pc2_1 = 0
        IFUCube.meta.wcsinfo.pc2_2 = 1
        IFUCube.meta.wcsinfo.pc2_3 = 0

        IFUCube.meta.wcsinfo.pc3_1 = 0
        IFUCube.meta.wcsinfo.pc3_2 = 0
        IFUCube.meta.wcsinfo.pc3_3 = 1

        IFUCube.meta.ifu.flux_extension = 'SCI'
        IFUCube.meta.ifu.error_extension = 'ERR'
        IFUCube.meta.ifu.error_type = 'ERR'
        IFUCube.meta.ifu.dq_extension = 'DQ'
        IFUCube.meta.ifu.roi_spatial = self.rois
        IFUCube.meta.ifu.roi_wave = self.roiw
        IFUCube.meta.ifu.weighting = self.weighting
        IFUCube.meta.ifu.weight_power = self.weight_power

        with datamodels.open(self.input_models[j]) as input:
            IFUCube.meta.bunit_data = input.meta.bunit_data
            IFUCube.meta.bunit_err = input.meta.bunit_err

        if self.coord_system == 'alpha-beta' :
            IFUCube.meta.wcsinfo.cunit1 = 'arcsec'
            IFUCube.meta.wcsinfo.cunit2 = 'arcsec'

# we only need to check list_par1[0] and list_par2[0] because these types
# of cubes are made from 1 exposures (setup_cube checks this at the start
# of cube_build).
            if self.list_par1[0] == '1' and self.list_par2[0] == 'SHORT' :
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL1A'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE1A'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'SHORT':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL2A'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE2A'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'SHORT':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL3A'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE3A'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'SHORT':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL4A'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE4A'

            if self.list_par1[0] == '1' and self.list_par2[0] == 'MEDIUM':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL1B'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE1B'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'MEDIUM':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL2B'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE2B'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'MEDIUM':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL3B'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE3B'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'MEDIUM':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL4B'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE4B'

            if self.list_par1[0] == '1' and self.list_par2[0] == 'LONG':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL1C'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE1C'
            if self.list_par1[0] == '2' and self.list_par2[0] == 'LONG':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL2C'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE2C'
            if self.list_par1[0] == '3' and self.list_par2[0] == 'LONG':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL3C'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE3C'
            if self.list_par1[0] == '4' and self.list_par2[0] == 'LONG':
                IFUCube.meta.wcsinfo.ctype1 = 'MRSAL4C'
                IFUCube.meta.wcsinfo.ctype2 = 'MRSBE4C'


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
#********************************************************************************
    def find_output_type(self):

        ValidChannel = ['1', '2', '3', '4']
        ValidSubchannel = ['SHORT', 'MEDIUM', 'LONG']

        nchannels = len(ValidChannel)
        nsubchannels = len(ValidSubchannel)


        ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H',
                    'G395M', 'G395H', 'PRISM']
        ValidFWA = ['F070LP', 'F070LP', 'F100LP', 'F100LP', 'F170LP',
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']

        nbands = len(ValidFWA)

#********************************************************************************
    def blend_output_metadata(self, IFUCube):

        """Create new output metadata based on blending all input metadata."""
        # Run fitsblender on output product
        output_file = IFUCube.meta.filename

#        log.info('Blending metadata for {}'.format(output_file))
        blendmeta.blendmodels(IFUCube, inputs=self.input_models,
                              output=output_file)
