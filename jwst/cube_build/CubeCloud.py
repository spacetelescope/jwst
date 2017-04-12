
 # Routines used in Spectral Cube Building
from __future__ import absolute_import, print_function

import sys
import numpy as np
import math
import logging
from .. import datamodels
from ..assign_wcs import nirspec
from ..datamodels import dqflags
from . import cube
from . import coord
from gwcs import wcstools
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
#________________________________________________________________________________

def MatchDet2Cube(self, input_model,
                       x, y, file_slice_no,
                       Cube,spaxel,
                       c1_offset, c2_offset):
    """
    Short Summary
    -------------
    map x,y to Point cloud  in final coordinate system (xi,eta of cube)

    Parameters
    ----------
    x,y list of x and y values to map
    input_model: slope image
    file_no: the index on the files that are used to construct the Cube
    Cube: holds the basic information on the Cube (including wcs of Cube)
    v2v32radec: temporary (until information is contained in assign_wcs)
                holds the information to do the transformation from v2-v3 to ra-dec
    c1_offset, c2_offset: dither offsets for each file (default = 0)
    provided by the user
    islice : a NIRSPEC parameter - slice number

    Returns
    -------
    spaxel class matched to detector pixels with flux and weighting updated for each
    match


    """

#________________________________________________________________________________
    if(Cube.instrument == 'MIRI'):
        file_no = file_slice_no
        det2ab_transform = input_model.meta.wcs.get_transform('detector','alpha_beta')
        detector2v23 = input_model.meta.wcs.get_transform('detector', 'v2v3')
        v23toworld = input_model.meta.wcs.get_transform("v2v3","world")
        worldtov23 = input_model.meta.wcs.get_transform("world","v2v3")
        v2ab_transform = input_model.meta.wcs.get_transform('v2v3',
                                                            'alpha_beta')

        alpha, beta, wave = det2ab_transform(x, y)
        v2, v3, lam23 = detector2v23(x, y)
        ra,dec,lam = v23toworld(v2,v3,lam23)

        valid1 = np.isfinite(v2)
        valid2 = np.isfinite(v3)

        if(self.weighting == 'miripsf'):
            a = Cube.a_wave[file_no]
            c = Cube.c_wave[file_no]
            wa = Cube.a_weight[file_no]
            wc = Cube.c_weight[file_no]

            # transform Cube Spaxel centers to alpha,beta of exposure
            # for MIRI weighting parameters are based on distance in alpha-beta coord system
            # transform the cube coordinate values to alpha and beta values
            # xi,eta -> ra,dec
            # ra-dec -> v2,v3
            # v2,v3 -> local alph,beta

    elif(Cube.instrument == 'NIRSPEC'):
        islice = file_slice_no
        slice_wcs = nirspec.nrs_wcs_set_input(input_model, islice)
        yrange = slice_wcs.bounding_box[1][0],slice_wcs.bounding_box[1][1]
        xrange = slice_wcs.bounding_box[0][0],slice_wcs.bounding_box[0][1]
        x,y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box, step=(1, 1), center=True)
        ra, dec, lam = slice_wcs(x, y) # return v2,v3 are in degrees

        valid1 = np.isfinite(ra)
        valid2 = np.isfinite(dec)
#________________________________________________________________________________
#________________________________________________________________________________
# Slices are curved on detector. A slice region is grabbed by corner regions so
# the region returned may include pixels not value for slice. There are gaps between
# the slices. Pixels not belonging to a slice are assigned NaN values.

    flux_all = input_model.data[y, x]
    error_all = input_model.err[y, x]
    dq_all = input_model.dq[y,x]

    valid3 = np.isfinite(lam)
    valid4 = np.isfinite(flux_all)
    valid = valid1 & valid2 & valid3 &valid4
#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded
# from the cube mapping
    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] +
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] +
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])

    # find the location of all the values to reject in cube building
    good_data = np.where((np.bitwise_and(dq_all, all_flags)==0) & (valid == True))

    # good data holds the location of pixels we want to map to cube
    flux = flux_all[good_data]
    error = error_all[good_data]
    wave = lam[good_data]

    xpix = x[good_data] # only used for testing
    ypix = y[good_data] # only used for testing

    ra = ra - c1_offset/3600.0
    dec = dec - c2_offset/3600.0
    ra_use = ra[good_data]
    dec_use = dec[good_data]
    if(Cube.instrument == 'MIRI'):
        # need alpha,beta if weigthing is miripsf or cubes in alpha-beta space
        alpha_det = alpha[good_data]
        beta_det = beta[good_data]
# MIRI can make cubes in alpha-beta:
    if(self.coord_system == 'alpha-beta'):
        coord1 = alpha[good_data]
        coord2 = beta[good_data]

    else:
        xi,eta = coord.radec2std(Cube.Crval1, Cube.Crval2,ra_use,dec_use) # xi,eta in arc seconds
        coord1 = xi
        coord2 = eta


#    print('ra use',ra_use[0:10])
#    print('dec use',dec_use[0:10])
#    print('ra use shape',ra_use.shape)
#        index_xy = np.nonzero( (xpix == 470) & (ypix ==1))[0]
#        print('wavelength',lam_use[index_xy])
#        mm = (lam_use[index_xy])
#        print('Wavelength at 470,1',mm)
#        n= len(lam)
#        for ii in range(0, n - 1):
#            if(lam[ii] < mm):
#                mnew = lam[ii]
#                print(mnew,mm,xpix[ii],ypix[ii],dq_all[ii])


    nxc = len(Cube.xcoord)
    nzc = len(Cube.zcoord)
    nyc = len(Cube.ycoord)
    nplane = Cube.naxis1 * Cube.naxis2
    lower_limit = 0.01

    iprint = 0
# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn  = coord1.size
#________________________________________________________________________________
    for ipt in range(0, nn - 1):

#________________________________________________________________________________
        # Coord1 and Coord2 are in the coordinate system of the cube.
        # using the Cube regularily spaced arrays - Cube.zcoord, xcoord,ycoord
        # find the spaxels that fall withing ROI of point cloud defined  by
        # coord1,coord2,wave

        indexz = np.where(abs(Cube.zcoord - wave[ipt]) <= self.roiw)
        indexx = np.where(abs(Cube.xcoord - coord1[ipt]) <= self.roi1)
        indexy = np.where(abs(Cube.ycoord - coord2[ipt]) <= self.roi2)

        zlam = Cube.zcoord[indexz] # Cube values for lambda axis
        xi = Cube.xcoord[indexx]   # Cube values for xi vector axis
        eta = Cube.ycoord[indexy]  # Cube values for eta vector axis

#        print('Cube Zcoord',Cube.zcoord[0:50])
#        print('wave',wave[ipt])
#        print(indexz[0])
#________________________________________________________________________________
# loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for iy, yy in enumerate(indexy[0]):
                for ix, xx in enumerate(indexx[0]):
#________________________________________________________________________________
                    if(self.weighting =='standard'):
                        d1 = (xi[ix] - coord1[ipt])/Cube.Cdelt1
                        d2 = (eta[iy] - coord2[ipt])/Cube.Cdelt2
                        d3 = (zlam[iz] - wave[ipt])/Cube.Cdelt3
                        weight_distance = math.sqrt(d1*d1 + d2*d2 + d3*d3)
                        weight_distance = math.pow(weight_distance,self.weight_power)
#________________________________________________________________________________
# if weight is miripsf -distances determined in alpha-beta coordinate system

                    elif(self.weighting =='miripsf'):

                        weights = FindNormalizationWeights(wave[ipt], a, c, wa, wc)
                        weight_alpha = weights[0]
                        weight_beta = weights[1]
                        weight_wave = weights[2]


                        ra_spaxel,dec_spaxel=coord.std2radec(Cube.Crval1,
                                                             Cube.Crval2,
                                                             xi[ix],eta[iy])


                        v2_spaxel,v3_spaxel,zl = worldtov23(ra_spaxel,dec_spaxel,zlam[iz])

                        alpha_spaxel,beta_spaxel,wave_spaxel = v2ab_transform(v2_spaxel,
                                                                              v3_spaxel,
                                                                              zlam[iz])
                        alpha_distance =alpha_det[ipt]-alpha_spaxel
                        beta_distance = beta_det[ipt]-beta_spaxel
                        wave_distance  = abs(wave[ipt]-wave_spaxel)

                        xn = alpha_distance/weight_alpha
                        yn = beta_distance/weight_beta
                        wn = wave_distance/weight_wave

                        # only included the spatial dimensions
                        weight_distance = math.sqrt(xn*xn + yn*yn + wn*wn)
                        weight_distance = math.pow(weight_distance,self.weight_power)
#________________________________________________________________________________
# We have found the weight_distance based on instrument type

                    if(weight_distance < lower_limit): weight_distance = lower_limit
                    weight_distance = 1.0 / weight_distance

                    cube_index = istart + yy * Cube.naxis1 + xx
                    spaxel[cube_index].flux = spaxel[cube_index].flux + weight_distance*flux[ipt]
                    spaxel[cube_index].flux_weight = spaxel[cube_index].flux_weight + weight_distance
                    spaxel[cube_index].iflux = spaxel[cube_index].iflux + 1


                    if(self.debug_pixel == 1 and self.xdebug == xx and
                       self.ydebug == yy and self.zdebug == zz ):

                        log.debug('For spaxel %d %d %d, detector x,y,flux %d %d %f %d %f '
                                  %(self.xdebug+1,self.ydebug+1,
                                    self.zdebug+1,xpix[ipt],ypix[ipt],
                                    flux[ipt],file_slice_no,weight_distance))
                        self.spaxel_debug.write('For spaxel %d %d %d, detector x,y,flux %d %d %f %d %f '
                                  %(self.xdebug+1,self.ydebug+1,
                                    self.zdebug+1,xpix[ipt],ypix[ipt],
                                    flux[ipt],file_slice_no,weight_distance) +' \n')
        iprint = iprint +1
        if(iprint == 10000):
            log.debug('Mapping point and finding ROI for point cloud # %d' %(ipt))
            iprint = 0
#        sys.exit('STOP')
#_______________________________________________________________________
def FindWaveWeights(channel, subchannel):
    """
    Short Summary
    -------------
    Get the wavelength normalization weights that we will use to normalize wavelengths.

    Parameters
    ----------
    channel- channel for point
    subchannel- subchannel for point

    Returns
    -------
    normalized weighting for wavelength for this channel, subchannel

    """

    if(channel == '1'):
        if(subchannel == 'SHORT'):
            a = 3050.0
            c = 3340.0
            wa = 4.91
            wc = 5.79
        elif(subchannel == 'MEDIUM'):
            a = 2920.0
            c = 3400.0
            wa = 5.6
            wc = 6.62
        elif(subchannel == 'LONG'):
            a = 2800.0
            c = 3220.0
            wa = 6.46
            wc = 7.63


    if(channel == '2'):
        if(subchannel == 'SHORT'):
            a = 2700.0
            c = 2800.0
            wa = 7.55
            wc = 8.91
        elif(subchannel == 'MEDIUM'):
            a = 2600.0
            c = 2880.0
            wa = 8.71
            wc = 10.34
        elif(subchannel == 'LONG'):
            a = 2590.0
            c = 3000.0
            wa = 9.89
            wc = 11.71

    if(channel == '3'):
        if(subchannel == 'SHORT'):
            a = 2390.0
            c = 2650.0
            wa = 11.50
            wc = 13.59
        elif(subchannel == 'MEDIUM'):
            a = 1600.0
            c = 2400.0
            wa = 13.19
            wc = 15.58
        elif(subchannel == 'LONG'):
            a = 1850.0
            c = 2550.0
            wa = 15.40
            wc = 18.14

    if(channel == '4'):
        if(subchannel == 'SHORT'):
            a = 1320.0
            c = 1720.0
            wa = 17.88
            wc = 21.34
        elif(subchannel == 'MEDIUM'):
            a = 1550.0
            c = 1600.0
            wa = 20.69
            wc = 24.68
        elif(subchannel == 'LONG'):
            a = 1450.0
            c = 1200.0
            wa = 23.83
            wc = 28.43
    return a, c, wa, wc

#_______________________________________________________________________

def FindNormalizationWeights(wavelength, a, c, wa, wc):
    """
    Short Summary
    -------------
    we need to normalize how to each point cloud in the spaxel. The normalization of
gs     weighting is determined from width of PSF as well as wavelength resolution

    Parameters
    ----------
    channel- channel for point
    subchannel- subchannel for point
    wavelength of point

    Returns
    -------
    normalized weighting for 3 dimension

    """
    alpha_weight = 1.0
    beta_weight = 1.0
    lambda_weight = 1.0

    beta_weight = 0.31 * (wavelength / 8.0)

    if(wavelength < 8.0):
        alpha_weight = 0.31
    else:
        alpha_weight = beta_weight


        # linear interpolation

    if (wavelength >= wa and wavelength <= wc):
        b = a + (c - a) * (wavelength - wa) / (wc - wa)

    elif (wavelength < wa):
        b = a
    else:
        b = c


    lambda_weight = wavelength / b

    weight = [alpha_weight, beta_weight, lambda_weight]
    return weight

