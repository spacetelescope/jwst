# Routines used for building cubes
from __future__ import absolute_import, print_function

import sys
import time
import numpy as np
import math
import json

from astropy.io import fits

from gwcs.utils import _domain_to_bounds
from ..associations import Association
from .. import datamodels
from ..assign_wcs import nirspec
from . import cube
from . import CubeOverlap
from . import CubeCloud
from . import coord

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


#********************************************************************************
def DetermineCubeCoverage(self, MasterTable):
#********************************************************************************
    """
    Short Summary
    -------------
    Function to determine which files contain channels and subchannels are used 
    in the creation of the cubes.
    For MIRI The channels  to be used are set by the association and the 
    subchannels are  determined from the data

    Parameter
    ----------
    self containing user set input parameters:
    self.channel, self.subchannel

    Returns
    -------
    fills in self.metadata['subchannel'] and self.metadata['channel']

    """
#________________________________________________________________________________
# loop over the file names

    if(self.metadata['instrument'] == 'MIRI'):
        ValidChannel = ['1', '2', '3', '4']
        ValidSubchannel = ['SHORT', 'MEDIUM', 'LONG']
        nchannels = len(ValidChannel)
        nsubchannels = len(ValidSubchannel)

        # for MIRI we can set the channel
        channellist = self.channel.split()
        user_clen = len(channellist)

        # the user has given the channel information
        if(user_clen != 0 and self.input_table_type == 'singleton'): 
            for j in range(user_clen):

                if channellist[j] in ValidChannel:
                    self.metadata['channel'].append(channellist[j])
                else:
                    log.error(' Invalid Channel give %s', channellist[j])

            self.metadata['channel'] = list(set(self.metadata['channel'])) # remove duplicates if needed

        for i in range(nchannels):
            for j in range(nsubchannels):
                nfiles = len(MasterTable.FileMap['MIRI'][ValidChannel[i]][ValidSubchannel[j]])
                if(nfiles > 0):
                    self.metadata['subchannel'].append(ValidSubchannel[j])
                    if(self.input_table_type == 'singleton' and user_clen == 0):
#usually filled in from reading assoication table
                        self.metadata['channel'].append(ValidChannel[i]) 

        self.metadata['subchannel'] = list(set(self.metadata['subchannel']))
        log.info('The desired cubes covers the MIRI Channels: %s', 
                 self.metadata['channel'])
        log.info('The desried cubes covers the MIRI subchannels: %s', 
                 self.metadata['subchannel'])


        number_channels = len(self.metadata['channel'])
        number_subchannels = len(self.metadata['subchannel'])

        if(number_channels == 0):
            raise ErrorNoChannels(
                "The cube  does not cover any channels, change parameter channel")
        if(number_subchannels == 0):
            raise ErrorNoSubchannels(
                "The cube  does not cover any subchannels, change parameter subchannel")
#______________________________________________________________________
    if(self.metadata['instrument'] == 'NIRSPEC'):

        ValidFWA = ['F070LP', 'F070LP', 'F100LP', 'F100LP', 'F170LP', 
                    'F170LP', 'F290LP', 'F290LP', 'CLEAR']
        ValidGWA = ['G140M', 'G140H', 'G140M', 'G140H', 'G235M', 'G235H', 
                    'G395M', 'G395H', 'PRISM']
        ntypes = len(ValidFWA)

        for j in range(ntypes):
            nfiles = len(MasterTable.FileMap['NIRSPEC'][ValidFWA[j]][ValidGWA[j]])
            if(nfiles > 0):
                self.metadata['filter'].append(ValidFWA[j])
                self.metadata['grating'].append(ValidGWA[j])

        self.metadata['filter'] = list(set(self.metadata['filter']))
        self.metadata['grating'] = list(set(self.metadata['grating']))
        log.debug('The desired cubes covers the NIRSPEC FWA  %s', 
                  self.metadata['filter'])
        log.debug('The desried cubes covers the NIRSPEC GWA: %s', 
                  self.metadata['grating'])

        number_filters = len(self.metadata['filter'])
        number_gratings = len(self.metadata['grating'])

        if(number_filters == 0):
            raise ErrorNoChannels("The cube  does not cover any filters")
        if(number_gratings == 0):
            raise ErrorNoSubchannels("The cube  does not cover any gratings")


#********************************************************************************
def FilesinCube(self, input_table, iproduct, MasterTable):
#********************************************************************************
    """
    Short Summary
    -------------
    Fill in the MasterTable which holds the files that the cube will be constructed 
    from. Since MIRI has 2 channels per image this MASTERTable helps to figure out
    which data needs to be use.
    THe MasterTable for MIRI is broken down by channel and subchannel.
    For each channel/subchannel combination - a file is listed that covers those options
    For NIRSPEC the table contains the Grating and Filter for each file. 

    If there is a dither offet file then the master table also holds the ra,dec offset for
    each file.

    Parameters
    ----------
    input_table: Association Table
    iproduct: index for the product in the assocation table. The product the is top grouping
    of how the data is organized. The Mastertable only holds the data for 1 product at a time.


    Returns
    -------
    MasterTable filled in with files needed
    num: number of files to create cube from
    detector

    """
    detector = ''
    input_filenames = []
    input_models = []
    i = 0
    num = 0

    if len(input_table.input_models) > 0:  # this is a single file

        input_models.append(input_table.input_models)
        input_filenames.append(input_table.filename)
        num = 1
    else:

        #channels = input_table.asn_table['products'][iproduct]['ch']
        channels = ['1']
        channellist = list(channels)
        num_ch = len(channellist)
        ValidChannel = ['1', '2', '3', '4']

        for j in range(num_ch):
            if channellist[j] in ValidChannel:
                self.metadata['channel'].append(channellist[j])
            else:
                log.error(' Invalid Channel %s', channellist[j])

        for m in input_table.asn_table['products'][iproduct]['members']:

            input_filenames.append(m['expname'])

            i = i + 1

    num = len(input_filenames)

#________________________________________________________________________________
# Loop over input list of files
    for i in range(num):

        ifile = input_filenames[i]

        # Open the input data model
        # Fill in the FileMap information

        with datamodels.ImageModel(ifile) as input_model:

            detector = input_model.meta.instrument.detector
            instrument = input_model.meta.instrument.name
            assign_wcs = input_model.meta.cal_step.assign_wcs

            if(assign_wcs != 'COMPLETE'):
                raise ErrorNoAssignWCS("Assign WCS has not been run on file %s", ifile)

            #________________________________________________________________________________
            #MIRI instrument
            #________________________________________________________________________________

            if instrument == 'MIRI':
                channel = input_model.meta.instrument.channel
                subchannel = input_model.meta.instrument.band
            #________________________________________________________________________________
                clenf = len(channel)
                for k in range(clenf):
                    MasterTable.FileMap['MIRI'][channel[k]][subchannel].append(ifile)
                    ioffset = len(self.ra_offset)
                    if (ioffset > 0):
                        ra_offset = self.ra_offset[i]
                        dec_offset = self.dec_offset[i]
                        MasterTable.FileOffset[channel[k]][subchannel]['C1'].append(ra_offset)
                        MasterTable.FileOffset[channel[k]][subchannel]['C2'].append(dec_offset)
            #________________________________________________________________________________
            elif instrument== 'NIRSPEC':

                fwa = input_model.meta.instrument.filter
                gwa = input_model.meta.instrument.grating

                MasterTable.FileMap['NIRSPEC'][fwa][gwa].append(ifile)
            else:

                print('Instrument not valid for cube')

    return num, instrument,detector
#********************************************************************************

def DetermineScale(Cube, InstrumentInfo):
#********************************************************************************
    """
    Short Summary
    -------------
    Determine the scale (sampling) in the 3 dimensions for the cube

    Parameters
    ----------
    Cube: Class holding basic information on cube
    InstrumentInfo holds the defaults scales for each channel/subchannel

    Returns
    -------
    scale, holding the scale for the 3 dimensions of the cube/

    """
    a = Cube.detector
    scale = [0, 0, 0]

    if(Cube.instrument == 'MIRI'):
        number_channels = len(Cube.channel)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_channels):
            this_channel = Cube.channel[i]
            a_scale, b_scale, wscale = InstrumentInfo.GetScale(this_channel)

            if(a_scale < min_a):
                min_a = a_scale
            if(b_scale < min_b):
                min_b = b_scale
            if(wscale < min_w):
                min_w = wscale

        scale = [min_a, min_b, min_w]

    elif(Cube.instrument == 'NIRSPEC'):
        number_gratings = len(Cube.grating)
        min_a = 1000.00
        min_b = 1000.00
        min_w = 1000.00

        for i in range(number_gratings):
            this_gwa = Cube.grating[i]
            a_scale, b_scale, wscale = InstrumentInfo.GetScale(this_gwa)

            if(a_scale < min_a):
                min_a = a_scale
            if(b_scale < min_b):
                min_b = b_scale
            if(wscale < min_w):
                min_w = wscale

        scale = [min_a, min_b, min_w]

    return scale
#_______________________________________________________________________


#********************************************************************************
def FindFootPrintMIRI(self, input, this_channel, InstrumentInfo):
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
    spaxial coordinates are in units of arc secons. 
    """
    # x,y values for channel - convert to output coordinate system
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers

    xstart, xend = InstrumentInfo.GetMIRISliceEndPts(this_channel)

    y, x = np.mgrid[:1024, xstart:xend]

    coord1 = np.zeros(y.shape)
    coord2 = np.zeros(y.shape)
    lam = np.zeros(y.shape)

    if (self.coord_system == 'alpha-beta'):
        detector2alpha_beta = input.meta.wcs.get_transform('detector', 'alpha_beta')
        coord1, coord2, lam = detector2alpha_beta(x, y)

    elif (self.coord_system == 'ra-dec'):
        detector2v23 = input.meta.wcs.get_transform('detector', 'V2_V3')
        #coord1,coord2,lam = input.meta.wcs(x,y)

        v2, v3, lam = detector2v23(x, y) # v2,v3 are in units of arc minutes 
        ra_ref = input.meta.wcsinfo.ra_ref # degrees
        dec_ref = input.meta.wcsinfo.dec_ref # degrees
        roll_ref = input.meta.wcsinfo.roll_ref # degrees 
        v2_ref = input.meta.wcsinfo.v2_ref # arc min
        v3_ref = input.meta.wcsinfo.v3_ref # arc min         
        coord1,coord2 = coord.V2V32RADEC(ra_ref,dec_ref,roll_ref,v2_ref,v3_ref,
                                         v2,v3) # return ra and dec in degrees

    else:
        # error the coordinate system is not defined
        raise NoCoordSystem(" The output cube coordinate system is not definded")

    a_min = np.nanmin(coord1)
    a_max = np.nanmax(coord1)

    b_min = np.nanmin(coord2)
    b_max = np.nanmax(coord2)

    lambda_min = np.nanmin(lam)
    lambda_max = np.nanmax(lam)

    print('return from footprint',a_min,a_max,b_min,b_max)

    
    return a_min, a_max, b_min, b_max, lambda_min, lambda_max


#********************************************************************************
def FindFootPrintNIRSPEC(self, input, this_channel):
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
    this_channel: channel working with


    Returns
    -------
    min and max spaxial coordinates  and wavelength for channel.

    """
    # loop over all the region (Slices) in the Channel
    # based on regions mask (indexed by slice number) find all the detector
    # x,y values for slice. Then convert the x,y values to  v2,v3,lambda
    # return the min & max of spatial coords and wavelength  - these are of the pixel centers



    start_slice = 0
    end_slice = 29

    nslices = end_slice - start_slice + 1

    a_slice = np.zeros(nslices * 2)
    b_slice = np.zeros(nslices * 2)
    lambda_slice = np.zeros(nslices * 2)

    regions = list(range(start_slice, end_slice + 1))
    k = 0



    # for NIRSPEC there are 30 regions
    for i in regions:

        slice_wcs = nirspec.nrs_wcs_set_input(input, 0, i)

        yrange = slice_wcs.domain[1]['lower'],slice_wcs.domain[1]['upper']
        xrange = slice_wcs.domain[0]['lower'],slice_wcs.domain[0]['upper']
        y, x = np.mgrid[yrange[0]:yrange[1], xrange[0]:xrange[1]]
        v2, v3, lam = slice_wcs(x, y) # return v2,v3 are in degrees

        print(yrange)
        print(xrange)
#        print('v2 ',v2.shape,v2)
#        print('v3 ',v3)


        v2 = v2 * 60.0 # convert to arc minutes
        v3 = v3 * 60.0 # convert to arc minutes


        ra_ref = 45.0 # degrees
        dec_ref = 0.0 # degrees
        roll_ref = 0.0 # degrees 
        v2_ref = 300.2961/60.0 # arc min
        v3_ref = -497.9/60.0 # arc min         
        coord1,coord2 = coord.V2V32RADEC(ra_ref,dec_ref,roll_ref,v2_ref,v3_ref,
                                         v2,v3) # return ra and dec in degrees


        a_slice[k] = np.nanmin(coord1)
        a_slice[k + 1] = np.nanmax(coord1)

        b_slice[k] = np.nanmin(coord2)
        b_slice[k + 1] = np.nanmax(coord2)

        lambda_slice[k] = np.nanmin(lam)
        lambda_slice[k + 1] = np.nanmax(lam)

        k = k + 2

    a_min = min(a_slice)
    a_max = max(a_slice)

    b_min = min(b_slice)
    b_max = max(b_slice)

    lambda_min = min(lambda_slice)
    lambda_max = max(lambda_slice)

    print('max a',a_min,a_max, (a_max-a_min)*60.0)
    print('max b',b_min,b_max, (b_max-b_min)*60.0)
    print('wave',lambda_min,lambda_max)
    sys.exit('STOP')
    return a_min, a_max, b_min, b_max, lambda_min, lambda_max

#_______________________________________________________________________
#********************************************************************************
def DetermineCubeSize(self, Cube, MasterTable, InstrumentInfo):
#********************************************************************************
    """
    Short Summary
    -------------
    Function to determine the min and max coordinates of the spectral cube,given channel & subchannel

    Parameter
    ----------
    Cube: class the holds the basic paramters of the IFU cube to be created
    MasterTable:  A table that contains the channel/subchannel or filter/grating for each input file
    InstrumentInfo: Default information on the MIRI and NIRSPEC instruments. This information might
                    contained in a different file in the future. Probably a reference file

    Returns
    -------
    Cube Dimension Information:

    Footprint of cube: min and max of coordinates of cube. If an offset list is provided then these values are applied.
    if the coordinate system is alpha-beta (MIRI) then min and max coordinates of alpha (arc sec),
    beta (arc sec) and lambda (microns) 
    if the coordinate system is ra-dec then the min and max of ra(degress), dec (degrees) and lambda (microns)
    is returned. 


    """
    instrument = Cube.instrument
    parameter1 = list()
    parameter2 = list()
    if(instrument == 'MIRI'):
        parameter1 = Cube.channel
        parameter2 = Cube.subchannel
    elif(instrument == 'NIRSPEC'):
        parameter1 = Cube.filter
        parameter2 = Cube.grating


    number_par1 = len(parameter1)
    number_par2 = len(parameter2)

    a_min = list()
    a_max = list()
    b_min = list()
    b_max = list()
    lambda_min = list()
    lambda_max = list()

    for i in range(number_par1):
        this_a = parameter1[i]
        for j in range(number_par2):
            this_b = parameter2[j]
            n = len(MasterTable.FileMap[instrument][this_a][this_b])
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

                ifile = MasterTable.FileMap[instrument][this_a][this_b][k]
                ioffset = len(MasterTable.FileOffset[this_a][this_b]['C1'])
                if(ioffset == n):
                    c1_offset = MasterTable.FileOffset[this_a][this_b]['C1'][k]
                    c2_offset = MasterTable.FileOffset[this_a][this_b]['C2'][k]
                    print('offset to apply in arcseonds (ra,dec) ', c1_offset, c2_offset)
#________________________________________________________________________________

                # Open the input data model
                with datamodels.ImageModel(ifile) as input_model:

                    t0 = time.time()
                    if(instrument == 'NIRSPEC'):
                        ChannelFootPrint = FindFootPrintNIRSPEC(self, input_model, this_a)
                        amin, amax, bmin, bmax, lmin, lmax = ChannelFootPrint
                        #print(amin,amax,bmin,bmax,lmin,lmax)
                        t1 = time.time()
                        log.info("Time find foot print = %.1f.s" % (t1 - t0,))
#________________________________________________________________________________
                    if(instrument == 'MIRI'):

                        ChannelFootPrint = FindFootPrintMIRI(self, input_model, this_a, InstrumentInfo)
                        t1 = time.time()
                        print("Time find foot print Quick = %.1f.s" % (t1 - t0,))
                        amin, amax, bmin, bmax, lmin, lmax = ChannelFootPrint


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

    CubeFootPrint = (final_a_min, final_a_max, final_b_min, final_b_max,
                     final_lambda_min, final_lambda_max)
    
    return CubeFootPrint
#________________________________________________________________________________


#********************************************************************************
def MapDetectorToCube(self, this_channel, this_subchannel, 
                      Cube, spaxel, 
                      PixelCloud, 
                      MasterTable, 
                      InstrumentInfo):
#********************************************************************************
    """
    Short Summary
    -------------
    Loop over files that cover the cube and map the detector pixel to Cube spaxels
    If dither offsets have been supplied then apply those values to the data

    Parameter
    ----------
    
    Cube - contains the basic header information of Cube
    spaxel: List of Spaxels

    Returns
    -------
    if(interpolation = area - only valid for alpha-beta
    or
    if(interpolation = pointcloud
    """

    instrument = Cube.instrument
    nfiles = len(MasterTable.FileMap[instrument][this_channel][this_subchannel])
    log.info('Number of files in cube %i', nfiles)

    # loop over the files that cover the spectral range the cube is for
    
    for k in range(nfiles):
        ifile = MasterTable.FileMap[instrument][this_channel][this_subchannel][k]
        #print(' On File k',k,nfiles)
        ioffset = len(MasterTable.FileOffset[this_channel][this_subchannel]['C1'])
        Cube.file.append(ifile)
        c1_offset = 0.0
        c2_offset = 0.0
        # c1_offset and c2_offset are the dither offset sets (in arc seconds)
        # by default these are zer0. The user has to supply these 
        if(ioffset == nfiles):
            c1_offset = MasterTable.FileOffset[this_channel][this_subchannel]['C1'][k]
            c2_offset = MasterTable.FileOffset[this_channel][this_subchannel]['C2'][k]

        log.info('Mapping file to cube %s ', ifile)

        # Open the input data model
        with datamodels.ImageModel(ifile) as input_model:
            det2ab_transform = input_model.meta.wcs.get_transform('detector', 'alpha_beta')
            v2ab_transform = input_model.meta.wcs.get_transform('V2_V3', 'alpha_beta')
            wave_weights = CubeCloud.FindWaveWeights(this_channel, this_subchannel)

            # for each file we need information that will be the same for all the pixels on the image
            # this information is used in the weight schemem on how to combine the surface brightness
            # information. The Cube class stores these paramters as a series of lists.  

            Cube.a_wave.append(wave_weights[0])
            Cube.c_wave.append(wave_weights[1])
            Cube.a_weight.append(wave_weights[2])
            Cube.c_weight.append(wave_weights[3])
            Cube.transform.append(v2ab_transform)

            # read in the V2-V3 to RA-DEC information
            ra_ref = input_model.meta.wcsinfo.ra_ref # degrees
            dec_ref = input_model.meta.wcsinfo.dec_ref # degrees
            roll_ref = input_model.meta.wcsinfo.roll_ref # degrees 
            v2_ref = input_model.meta.wcsinfo.v2_ref # arc min
            v3_ref = input_model.meta.wcsinfo.v3_ref # arc min 

            v2v32radec = ra_ref,dec_ref,roll_ref,v2_ref,v3_ref  # temporarily
            # store the info needed to transform to ra,dec (this will later
            # be in assign_wcs) 

            Cube.ra_ref.append(ra_ref)
            Cube.dec_ref.append(dec_ref)
            Cube.roll_ref.append(roll_ref)
            Cube.v2_ref.append(v2_ref)
            Cube.v3_ref.append(v3_ref)

#________________________________________________________________________________
# Standard method 
            if(self.interpolation == 'pointcloud'):
                xstart, xend = InstrumentInfo.GetMIRISliceEndPts(this_channel)
                y, x = np.mgrid[:1024, xstart:xend]
                y = np.reshape(y, y.size)
                x = np.reshape(x, x.size)

                t0 = time.time()
                cloud = CubeCloud.MakePointCloudMIRI(self,input_model,
                                                     x, y, k, 
                                                     Cube,
                                                     v2v32radec,
                                                     c1_offset, c2_offset)
                n = PixelCloud.size
                print('size of PixelCloud',n)
                if(n == 10):  # If first time
                    PixelCloud = cloud
                    print(' 1 pixelcloud',PixelCloud.shape, cloud.shape)
                else:    #  add information for another slice  to the  PixelCloud

                    PixelCloud = np.hstack((PixelCloud, cloud))
                print('2 pixelcloud',PixelCloud.shape, cloud.shape)
                t1 = time.time()
                log.debug("Time Map one Channel  to Cloud = %.1f.s" % (t1 - t0,))
#________________________________________________________________________________

#2D area method - only works for single files and coord_system = 'alpha-beta'
            if(self.interpolation == 'area'):

                start_region = InstrumentInfo.GetStartSlice(this_channel)
                end_region = InstrumentInfo.GetEndSlice(this_channel)
                regions = list(range(start_region, end_region + 1))

                for i in regions:
                    log.info('Working on Slice # %d', i)

                    y, x = (det2ab_transform.label_mapper.mapper == i).nonzero()

                    # spaxel object holds all needed information in a set of lists
                    #    flux (of overlapping detector pixel)
                    #    error (of overlapping detector pixel)
                           #    overlap ratio
                    #    beta distance

                    index = np.where(y < 1023) # getting pixel corner - ytop = y + 1 (routine fails for y = 1024)
                    y = y[index]
                    x = x[index]
                    t0 = time.time()

                    beta_width = Cube.Cdelt2
                    CubeOverlap.SpaxelOverlap(self, x, y, i, 
                                              start_region, 
                                              input_model, 
                                              det2ab_transform, 
                                              beta_width, 
                                              Cube, spaxel)
                    t1 = time.time()
                    log.debug("Time Map one Slice  to Cube = %.1f.s" % (t1 - t0,))

    print('pixelcloud',PixelCloud.shape)
    return PixelCloud

#********************************************************************************
def FindCubeFlux(self, Cube, spaxel, PixelCloud):
#********************************************************************************
    """
    Short Summary
    -------------
    Depending on the interpolation method, find the flux for each spaxel value

    Parameter
    ----------
    Cube - contains the basic header information of Cube
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
            s = len(spaxel[i].pixel_overlap)
            if(s > 0):
                CubeOverlap.SpaxelFlux(self.roi2, i, Cube, spaxel)

    elif self.interpolation == 'pointcloud':
        icube = 0
        t0 = time.time()
        iz = 0

        for z in Cube.zcoord:
            iy = 0

            for y in Cube.ycoord:
                ix = 0
                for x in Cube.xcoord:
                    num = len(spaxel[icube].ipointcloud)
                    if(num > 0):
                        pointcloud_index = spaxel[icube].ipointcloud
                        weightpt = spaxel[icube].pointcloud_weight
                        pixelflux = PixelCloud[5, pointcloud_index]

                        weight = 0
                        value = 0
                        for j in range(num):
                            #weightpt[j] = 1
                            weight = weight + weightpt[j]
                            value = value + weightpt[j] * pixelflux[j]

                            if(icube == 253984  ):
                                print('Checking ', icube, ix, iy, iz)
                                print('icube', icube)
                                print('pointcloud', pointcloud_index[j])
                                print('flux = {0:.5f}'.format(pixelflux[j]))
                                print('w', weightpt[j])
                                print(' ',weightpt[j] * pixelflux[j])
                                print('num',num)
                                
                        if(weight != 0):
                            value = value / weight
                            spaxel[icube].flux = value
                            if(icube == 253984  ):
                                print('Final Flux', value * weight, weight, value,num)


                    icube = icube + 1
                    ix = ix + 1
                iy = iy + 1
            iz = iz + 1

        t1 = time.time()
        log.info("Time to interpolate at spaxel values = %.1f.s" % (t1 - t0,))


#________________________________________________________________________________
# We might need a dither offset file. A file for each input image that defines the 
# ra and dec offset (in arc seconds) that the files are offset from each other. 
# For testing this is useful but possibily this might be useful during flight if the
# images need an additional offset applied to them 
def ReadOffSetFile(self):
    print('Going to read offset list', self.offset_list)
    f = open(self.offset_list, 'r')
    i = 0
    for line in f:
        #print('input offset', line)
        offset_str = line.split()
        offset = [float(xy) for xy in offset_str]
        ra_off = offset[0]
        dec_off = offset[1]

        self.ra_offset.append(ra_off)
        self.dec_offset.append(dec_off)
        print('Offset in V2,V3 (in arc seconds) for exposure', ra_off, dec_off, i)
        i = i + 1
    f.close()
#********************************************************************************
def UpdateOutPutName(self):

    ch_name = self.channel.replace(" ", "")
    newname = self.output_name + '-CH' + ch_name + '_IFUCube.fits'
    return newname


#********************************************************************************
def CheckCubeType(self):

    print(self.interpolation, len(self.metadata['channel']), self.metadata['channel'])
    if(self.interpolation == "area"):
        if(self.metadata['number_files'] > 1):
            raise IncorrectInput("For interpolation = area, only one file can be used to created the cube")

        if(len(self.metadata['channel']) > 1):
            raise IncorrectInput("For interpolation = area, only channel can be used to created the cube")

    if(self.coord_system == "alpha-beta"):
        if(self.metadata['number_files'] > 1):
            raise IncorrectInput("Cubes built in alpha-beta coordinate system are built from a single file")


#********************************************************************************
def WriteCube(self, Cube, spaxel):

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
    no return = writes file

    """
    #pull out data into array

    data = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    idata = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    dq_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))
    err_cube = np.zeros((Cube.naxis3, Cube.naxis2, Cube.naxis1))

    icube = 0
    for z in range(Cube.naxis3):
        for y in range(Cube.naxis2):
            for x in range(Cube.naxis1):
                data[z, y, x] = spaxel[icube].flux
                idata[z, y, x] = len(spaxel[icube].ipointcloud)

                icube = icube + 1
    name = Cube.output_name
    new_model = datamodels.IFUCubeModel(data=data, dq=dq_cube, err=err_cube, weightmap=idata)

    new_model.meta.wcsinfo.crval1 = Cube.Crval1
    new_model.meta.wcsinfo.crval2 = Cube.Crval2
    new_model.meta.wcsinfo.crval3 = Cube.Crval3
    new_model.meta.wcsinfo.crpix1 = Cube.Crpix1
    new_model.meta.wcsinfo.crpix2 = Cube.Crpix2
    new_model.meta.wcsinfo.crpix3 = Cube.Crpix3

    new_model.meta.wcsinfo.crdelt1 = Cube.Cdelt1
    new_model.meta.wcsinfo.crdelt2 = Cube.Cdelt2
    new_model.meta.wcsinfo.crdelt3 = Cube.Cdelt3

    new_model.meta.wcsinfo.ctype1 = 'RA---TAN'
    new_model.meta.wcsinfo.ctype2 = 'DEC--TAN'
    new_model.meta.wcsinfo.ctype3 = 'WAVE'

    new_model.meta.wcsinfo.cunit1 = 'DEG'
    new_model.meta.wcsinfo.cunit2 = 'DEG'
    new_model.meta.wcsinfo.cunit3 = 'MICRON'

#    new_model.meta.wcsinfo.waverange_start = 
#    new_model.meta.wcsinfo.waverange_end = 
    new_model.meta.flux_extension = 'SCI'
    new_model.meta.error_extension = 'ERR'
    new_model.meta.dq_extension = 'DQ'
    new_model.meta.weightmap = 'WMAP'
    new_model.error_type = 'ERR'


    new_model.save(name)
    new_model.close()


    log.info('Wrote %s', name)
    #return new_model

#********************************************************************************
class ErrorNoData(Exception):
    pass

class ErrorNoAssignWCS(Exception):
    pass

class ErrorNoChannels(Exception):
    pass

class ErrorNoSubchannels(Exception):
    pass

class ErrorNoCubes(Exception):
    pass

class IncorrectInput(Exception):
    pass

class NoCoordSystem(Exception):
    pass
#********************************************************************************
class IFUCubeInput(object):
#********************************************************************************

    """
    Class to handle reading the input to the processing, which
    can be a single science exposure or an IFU cube association table.
    The input and output member info is loaded into an ASN table model.
    """

    template = {"asn_rule": "",
              "target": "",
              "asn_pool": "",
              "asn_type": "",
              "products": [
                  {"name": "", "ch": "",
                   "members": [
                      {"exptype": "",
                       "expname": ""}
                      ]
                  }
                ]
              }

    def __init__(self, input):

        self.input = input # keep a record of original input name for later

        self.input_models = []

        print('input trying to open',input)
        if isinstance(input, datamodels.ImageModel):
            log.debug( 'going to read imagemodel')
            # It's a single image that's been passed in as a model
            self.interpret_image_model(input)
        elif isinstance(input, str):
            try:
                # The name of an association table
                with open(input, 'r') as input_fh:
                    self.asn_table = Association.load(input_fh)
            except:
                # The name of a single image file
                log.debug( 'going to read a single file')
                self.filename = input  # temp until figure out model.meta.filename
                self.interpret_image_model(datamodels.ImageModel(input))
            #except:
            #    raise IOError("Input is not single filei, fits file or association table")


    def interpret_image_model(self, model):
        """ Interpret image model as single member association data product.
        """

        # An in-memory ImageModel for a single exposure was provided as input
        self.asn_table = self.template
        self.asn_table['target'] = model.meta.target.catalog_name
        self.asn_table['asn_rule'] = 'singleton'
        self.asn_table['asn_type'] = 'singleton'

        self.asn_table['products'][0]['name'] = self.build_product_name(self.filename)
        self.rootname = self.filename[:self.filename.rfind('_')]
        self.asn_table['products'][0]['members'][0]['expname'] = self.filename
        self.input_models.append(model)

    def build_product_name(self, filename):
        indx = filename.rfind('.fits')
        single_product = filename[:indx]
        return single_product

    def get_inputs(self, product=0):
        members = []
        for p in self.asn_table['products'][product]['members']:
            members.append(p['expname'])
        return members
    def get_outputs(self, product=0):
        return self.asn_table['products'][product]['name']
