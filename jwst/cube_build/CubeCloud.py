# Routines used in Spectral Cube Building
from __future__ import absolute_import, print_function

import sys
import numpy as np
import math
from .. import datamodels
from ..assign_wcs import nirspec
from ..datamodels import dqflags
from . import cube
from . import coord
#________________________________________________________________________________

def MakePointCloudMIRI(self, input_model,
                       x, y, file_no, 
                       Cube,
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




    Returns
    -------
    coordinates of x,y in point cloud


    """
#________________________________________________________________________________
    det2ab_transform = input_model.meta.wcs.get_transform('detector','alpha_beta')
    detector2v23 = input_model.meta.wcs.get_transform('detector', 'v2v3')
    v23toworld = input_model.meta.wcs.get_transform("v2v3","world")
    worldtov23 = input_model.meta.wcs.get_transform("world","v2v3")

    alpha, beta, wave = det2ab_transform(x, y)
    v2, v3, lam = detector2v23(x, y)
    flux_all = input_model.data[y, x]
    error_all = input_model.err[y, x]
    dq_all = input_model.dq[y,x]

#________________________________________________________________________________
# in the slice gaps the v2,v3 and lam are NaN values. Select the valid slice pixels
# based on if v2,v3,lam are finite 
    valid1 = np.isfinite(v2) 
    valid2 = np.isfinite(v3)
    valid3 = np.isfinite(lam)  

    valid = dq_all.copy() * 0 
    index  = np.where(np.logical_and(np.logical_and(valid1,valid2),valid3))
    valid[index] = 1

#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded 
# from the cube mapping    
    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] + 
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] + 
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])

    # find the location of all the values to reject in cube building    
    good_data = np.where((np.bitwise_and(dq_all, all_flags)==0) & (valid == 1))

    # good data holds the location of pixels we want to map to cube 
    flux = flux_all[good_data]
    error = error_all[good_data]
    alpha = alpha[good_data]
    beta = beta[good_data]
    xpix = x[good_data] # only used for testing
    ypix = y[good_data] # only used for testing

#    ia = 0 
#    for f in flux:
#        if(f < -10):
#            print('very low',flux[ia],xpix[ia],ypix[ia])
#        ia = ia + 1
    if(self.coord_system == 'alpha-beta'):
        coord1 = alpha
        coord2 = beta
    else:
        v2_use = v2[good_data] #arc mins
        v3_use = v3[good_data] #arc mins
        lam_use = lam[good_data]

        ra,dec,wave = v23toworld(v2_use,v3_use,lam_use)

#        print('ra     ',ra[15:20])
#        print('dec     ',dec[15:20])

        ra = ra - c1_offset/3600.0
        dec = dec - c2_offset/3600.0
        xi,eta = coord.radec2std(Cube.Crval1, Cube.Crval2,ra,dec) # xi,eta in arc seconds
        coord1 = xi
        coord2 = eta

    ifile = np.zeros(flux.shape, dtype='int') + int(file_no)

    # get in form of 8 columns of data - shove the information in an array.
    # xpix,ypix used for testing
    cloud = np.asarray([coord1, coord2, wave, alpha, beta, flux, error, ifile, xpix, ypix])
#    cloud = np.asarray([coord1, coord2, wave, alpha, beta, flux, error, ifile])

    return cloud

#________________________________________________________________________________

def MakePointCloudNIRSPEC(self, input_model,
                          file_no,
                          islice,
                          Cube,
                          c1_offset, c2_offset):
    """

    Short Summary
    -------------
    For NIRSPEC IFU map x,y to Point cloud  in final coordinate system (xi,eta of cube) 

    Parameters
    ----------
    input_model: slope image
    file_no: the index on the files that are used to construct the Cube
    Cube: holds the basic information on the Cube (including wcs of Cube)
    v2v32radec: temporary (until information is contained in assign_wcs) 
                holds the information to do the transformation from v2-v3 to ra-dec
    c1_offset, c2_offset: dither offsets for each file (default = 0)
               provided by the user

    Returns
    -------
    each valid detector mapped to  in point cloud


    """
#________________________________________________________________________________

    slice_wcs = nirspec.nrs_wcs_set_input(input_model, islice)
    yrange = slice_wcs.domain[1]['lower'],slice_wcs.domain[1]['upper']
    xrange = slice_wcs.domain[0]['lower'],slice_wcs.domain[0]['upper']
    y, x = np.mgrid[yrange[0]:yrange[1], xrange[0]:xrange[1]]
    ra, dec, lam = slice_wcs(x, y) # return v2,v3 are in degrees

    print('yrange for slice',yrange,islice)
    print('xrange for slice',xrange,islice)

    flux_all = input_model.data[y, x]
    error_all = input_model.err[y, x]
    dq_all = input_model.dq[y,x]

#________________________________________________________________________________
# Slices are curved on detector. A slice region is grabbed by corner regions so
# the region returned may include pixels not value for slice 
    valid1 = np.isfinite(ra) 
    valid2 = np.isfinite(dec)
    valid3 = np.isfinite(lam)  

    valid = dq_all.copy() * 0 
    index  = np.where(np.logical_and(np.logical_and(valid1,valid2),valid3))
    valid[index] = 1

#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded 
# from the cube mapping    
    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] + 
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] + 
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])

    # find the location of all the values to reject in cube building    
    good_data = np.where((np.bitwise_and(dq_all, all_flags)==0) & (valid == 1))

    # good data holds the location of pixels we want to map to cube 
    flux = flux_all[good_data]
    error = error_all[good_data]
    alpha = flux*0
    beta = flux*0
    xpix = x[good_data] # only used for testing
    ypix = y[good_data] # only used for testing

    ra_use = ra[good_data] #arc mins
    dec_use = dec[good_data] #arc mins
    wave = lam[good_data]

    #print('shape of ra',len(ra),len(ra_use))


    ra_use = ra_use - c1_offset/3600.0
    dec_use = dec_use - c2_offset/3600.0
    xi,eta = coord.radec2std(Cube.Crval1, Cube.Crval2,ra_use,dec_use) # xi,eta in arc seconds
    coord1 = xi
    coord2 = eta

    ifile = np.zeros(flux.shape, dtype='int') + int(file_no)

    # stuff the point cloud arrays for this configuration into cloud 
    # Point cloud will eventually contain all the cloud values
    # xpix,ypix used for testing
    cloud = np.asarray([coord1, coord2, wave, alpha, beta, flux, error, ifile, xpix, ypix])

    return cloud
#______________________________________________________________________

def FindROI(self, Cube, spaxel, PointCloud):

    """
    Short Summary
    -------------
    using the point cloud loop over the Spaxel and find the point cloud members that
    fall withing the ROI of the spaxel center. 


    For MIRI the weighting of the Cloud points is based on the distance in the local
    MRS alpha-beta plane. Each cloud point as an associated alpha-beta coordinate
    The spaxel centers have xi,eta & V2,V3 so we need to know the channel and band 
    information and transform the V2,V3 coordinates back to alpha-beta

    Parameters
    ----------
    Cube: holds basic Cube information
    spaxel: a class that holds information on each spaxel in the cube. 
    PointCloud: array of point cloud members

    Returns
    -------
    location of x,y in Point Cloud as well as mapping of spaxel to each overlapping 
    PointCloud member


    """
#________________________________________________________________________________
    nxc = len(Cube.xcoord)
    nzc = len(Cube.zcoord)
    nyc = len(Cube.ycoord)

    nplane = Cube.naxis1 * Cube.naxis2
    lower_limit = 0.0001

    iprint = 0
    nn = len(PointCloud[0])
#    nn = 100

    print('number of elements in PT',nn)
    
# loop over each point cloud member - might want to change this to looping
# over spaxels but for now just keep it point cloud elements because it
# is easy to find ROI members because the cube spaxel values are regularily spaced
# and we can search over the vector of each axis of the cube rather than entire 
# point could
#________________________________________________________________________________
    for ipt in range(0, nn - 1):

        coord1 = PointCloud[0, ipt]  # Point cloud xi 
        coord2 = PointCloud[1, ipt]  # Point cloud etax
        wave = PointCloud[2,ipt]     # Point cloud wavelength 

        if(Cube.instrument == 'MIRI'):
            if(self.weighting == 'miripsf'):
                alpha = PointCloud[3, ipt]   
                beta = PointCloud[4, ipt]    

                ifile = int(PointCloud[7, ipt])
                a = Cube.a_wave[ifile]
                c = Cube.c_wave[ifile]
                wa = Cube.a_weight[ifile]
                wc = Cube.c_weight[ifile]
                weights = FindNormalizationWeights(wave, a, c, wa, wc)
                weight_alpha = weights[0]
                weight_beta = weights[1]
                weight_wave = weights[2]

        # transform Cube Spaxel centers to alpha,beta system
        # of point cloud member (only do this transformation for MIRI) 
        # for MIRI weighting parameters are based on distance in alpha-beta coord system
        # transform the cube coordinate values to alpha and beta values 
        # xi,eta -> ra,dec
        # ra-dec -> v2,v3 
        # v2,v3 -> local alph,beta

                v2ab_transform = Cube.transform_v23toab[ifile]
                worldtov23 = Cube.transform_worldtov23[ifile]

#________________________________________________________________________________

            if(self.coord_system == 'alpha-beta'):
                coord1 = alpha
                coord2 = beta
#________________________________________________________________________________        
        # Coord1 and Coord2 are in the coordinate system of the cube.
        # using the Cube regularily spaced arrays - Cube.zcoord, xcoord,ycoord
        # find the spaxels that fall withing ROI of point cloud

        # find values within the ROI

        indexz = np.where(abs(Cube.zcoord - wave) <= self.roiw)
        indexx = np.where(abs(Cube.xcoord - coord1) <= self.roi1)
        indexy = np.where(abs(Cube.ycoord - coord2) <= self.roi2)

        zlam = Cube.zcoord[indexz]
        xi = Cube.xcoord[indexx]   # Cube values for xi vector axis 
        eta = Cube.ycoord[indexy]  # Cube values for eta vector axis

#________________________________________________________________________________
# loop over the points in the ROI 
        iz = 0
        for zz in indexz[0]:
            istart = zz * nplane
            iy = 0
            for yy in indexy[0]:
                ix = 0
                for xx in indexx[0]:
#________________________________________________________________________________
# NIRSPEC instrument
                    # for NIRSPEC find distance between PT and Spaxel Center 
                    # in xi,eta coordinate system
                    if(Cube.instrument == 'NIRSPEC'):
                        d1 = abs(xi[ix] - coord1)
                        d2 = abs(eta[iy] - coord2)
                        weight_distance = d1*d1 + d2*d2 

#________________________________________________________________________________
# MIRI instrument
                    elif(Cube.instrument == 'MIRI'):
                        # weighting - standard - distance based on xi,eta distance
                        if(self.weighting =='standard'):
                            d1 = abs(xi[ix] - coord1)
                            d2 = abs(eta[iy] - coord2)
                            weight_distance = d1*d1 + d2*d2 

                        # For MIRI the distance between PT and Spaxel Center is
                        # in the alpha - beta cooridate system
                        elif(self.weighting =='miripsf'):

                            ra_spaxel,dec_spaxel=coord.std2radec(Cube.Crval1,
                                                             Cube.Crval2,
                                                             xi[ix],eta[iy])

 
                            v2_spaxel,v3_spaxel,zl = worldtov23(ra_spaxel,dec_spaxel,zlam[iz])
#                            print('v2      ',v2_spaxel*3600.0)
#                            print('v3      ',v3_spaxel*3600.0)

                            alpha_spaxel,beta_spaxel,wave_spaxel = v2ab_transform(v2_spaxel,
                                                                              v3_spaxel,
                                                                              zlam[iz])                    
#                            print('alpha_spaxel',alpha_spaxel)
#                            print('alpha       ',alpha)

#                            print('beta_spaxel ',beta_spaxel)
#                            print('beta        ',beta)
#                            sys.exit('STOP')
                            alpha_distance = abs(alpha-alpha_spaxel)
                            beta_distance = abs(beta-beta_spaxel)
                            wave_distance  = abs(wave-wave_spaxel)

                            xn = alpha_distance/weight_alpha
                            yn = beta_distance/weight_beta
                            wn = wave_distance/weight_wave
                                                         
                        # only included the spatial dimensions
                            weight_distance = xn*xn + yn*yn  
                                                          
#________________________________________________________________________________
# We have found the weight_distance based on instrument type

                    if(weight_distance < lower_limit): weight_distance = lower_limit
                    weight_distance = 1.0 / weight_distance

                    cube_index = istart + yy * Cube.naxis1 + xx
                    spaxel[cube_index].ipointcloud.append(ipt)
                    spaxel[cube_index].pointcloud_weight.append(weight_distance)

                    ix = ix + 1
#                    iprint = iprint + 1
#                    if(iprint > 80000):
#                        iprint = 0
#                        print('on point',ipt,nn)
                        
                
                iy = iy + 1
            iz = iz + 1


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



def FindNormalizationWeights(a, c, wa, wc, wavelength):
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

