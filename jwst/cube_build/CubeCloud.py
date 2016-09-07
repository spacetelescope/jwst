# Routines used in Spectral Cube Building
from __future__ import absolute_import, print_function

import sys
import numpy as np
import math
from .. import datamodels
from ..datamodels import dqflags
from . import CubeD2C
from . import cube
#________________________________________________________________________________


def MakePointCloudMIRI(self, x, y, file_no, c1_offset, c2_offset, input_model):
    """
    Short Summary
    -------------
    map x,y to Point cloud

    Parameters
    ----------
    x,y list of x and y values to map
    input_model: slope image
    transform: wsc transform to go from x,y to point cloud


    Returns
    -------
    location of x,y in Point Cloud as well as mapping of spaxel to each overlapping PointCloud member


    """

    det2ab_transform = input_model.meta.wcs.get_transform('detector', 'alpha_beta')
    alpha, beta, wave = det2ab_transform(x, y)

    detector2v23 = input_model.meta.wcs.get_transform('detector', 'V2_V3')

    v2, v3, lam = detector2v23(x, y)
    flux_all = input_model.data[y, x]
    error_all = input_model.err[y, x]
    dq_all = input_model.dq[y,x]
    #v2,v3,lam = input_model.meta.wcs(x,y)

    valid1 = np.isfinite(v2) 
    valid2 = np.isfinite(v3)
    valid3 = np.isfinite(lam)  
    #index = np.asarray(np.where(np.logical_and(np.logical_and(valid1,valid2),valid3)))


    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] + dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] + dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])
    
    good_data = np.asarray(np.where((np.bitwise_and(dq_all, all_flags)==0) & (np.logical_and(np.logical_and(valid1,valid2),valid3))))


    # find the location of all the values to reject in cube building
    #loc_do_not_use = np.bitwise_and(dq_all,dqflags.pixel['DO_NOT_USE'])  

    flux = flux_all[good_data]
    error = error_all[good_data]
    alpha = alpha[good_data]
    beta = beta[good_data]
    xpix = x[good_data]
    ypix = y[good_data]

    coord1 = v2[good_data] * 60.0
    coord2 = v3[good_data] * 60.0
    coord1 = coord1 - c1_offset 
    coord2 = coord2 - c2_offset 
    wave = lam[good_data]

    ifile = np.zeros(flux.shape, dtype='int') + int(file_no)

    # get in form of 8 columns of data - shove the information in an array.

    cloud = np.asarray([coord1, coord2, wave, alpha, beta, flux, error, ifile, xpix, ypix])

    return cloud
#______________________________________________________________________
def MakePointCloudMIRI_DistortionFile(self, x, y, file_no, c1_offset, c2_offset, channel, subchannel, sliceno, start_sliceno, input_model):


    """
    Short Summary
    -------------
    map x,y to Point cloud
    x,y detector values for sliceno

    Parameters
    ----------
    x,y list of x and y values to map
    channel: channel number working with
    subchannel: subchannel working with
    sliceno
    input_model: slope image
    transform: wsc transform to go from x,y to point cloud

    Returns
    -------
    location of x,y in Point Cloud as well as mapping of spaxel to each overlapping PointCloud member


    """
#________________________________________________________________________________
# if using distortion polynomials

    nn = len(x)
    sliceno_use = sliceno - start_sliceno + 1
    coord1 = list()
    coord2 = list()
    wave = list()
    alpha = list()
    beta = list()
    ifile = list()
    flux = list()
    error = list()
    x_pixel = list()
    y_pixel = list()
#________________________________________________________________________________
# loop over pixels in slice
#________________________________________________________________________________

    for ipixel in range(0, nn - 1):
        valid_pixel = True
        if(y[ipixel] >= 1024):
            valid_pixel = False
            #print(' Error ypixel = ',y[ipixel])

        if(valid_pixel):
            alpha_pixel, beta_pixel, wave_pixel = CubeD2C.xy2abl(self, sliceno_use - 1, x[ipixel], y[ipixel])

            flux_pixel = input_model.data[y[ipixel], x[ipixel]]
            error_pixel = input_model.err[y[ipixel], x[ipixel]]
            xan, yan = CubeD2C.ab2xyan(self, alpha_pixel, beta_pixel)
            v2, v3 = CubeD2C.xyan2v23(self, xan, yan)

            coord1_pixel = v2 * 60.0
            coord2_pixel = v3 * 60.0

            coord1_pixel = coord1_pixel + c1_offset / 60.0
            coord2_pixel = coord2_pixel + c2_offset / 60.0

            coord1.append(coord1_pixel)
            coord2.append(coord2_pixel)
            alpha.append(alpha_pixel)
            beta.append(beta_pixel)
            flux.append(flux_pixel)
            error.append(error_pixel)
            ifile.append(file_no)

            x_pixel.append(x[ipixel])
            y_pixel.append(y[ipixel])
    # get in form of 8 columns of data - shove the information in an array.
    #print('size',len(coord1),len(coord2))

    cloud = np.asarray([coord1, coord2, wave, alpha, beta, flux, error, ifile, x_pixel, y_pixel])
    return cloud

#_______________________________________________________________________



def FindROI(self, Cube, spaxel, PointCloud):

    nxc = len(Cube.xcoord)
    nzc = len(Cube.zcoord)
    nyc = len(Cube.ycoord)

    nplane = Cube.naxis1 * Cube.naxis2
    lower_limit = 0.0001

    iprint = 0
    nn = len(PointCloud[0])
#    print('number of elements in PT',nn)
    for ipt in range(0, nn - 1):

#        if(iprint == 0):
#            print('On pt member',ipt,'out of ',nn)
        ifile = int(PointCloud[7, ipt])

        a = Cube.a_wave[ifile]
        c = Cube.c_wave[ifile]
        wa = Cube.a_weight[ifile]
        wc = Cube.c_weight[ifile]
        wave = PointCloud[2, ipt]
        weights = FindNormalizationWeights(wave, a, c, wa, wc)

        weight_alpha = weights[0]
        weight_beta = weights[1]
        weight_wave = weights[2]

        coord1 = PointCloud[0, ipt]
        coord2 = PointCloud[1, ipt]
        alpha = PointCloud[3, ipt]
        beta = PointCloud[4, ipt]
        x = PointCloud[8, ipt]
        y = PointCloud[9, ipt]

        if(self.coord_system == 'alpha-beta'):
            coord1 = alpha
            coord2 = beta
# Map this point cloud to spaxel (distance from spaxel center = ROI)
# use the vector of cube centers in each dimension to speed things up
        # both coordinates are in arc seconds
        indexz = np.asarray(np.where(abs(Cube.zcoord - wave) <= self.roiw))
        indexx = np.asarray(np.where(abs(Cube.xcoord - coord1) <= self.roi1))
        indexy = np.asarray(np.where(abs(Cube.ycoord - coord2) <= self.roi2))


        # transform Cube Spaxel to alpha,beta system of Point Cloud
        # use inverse transfrom of V2,V3 back to local alpha,beta, lam
        # convert v2,v3 from arc seconds to arc min first


        zloc = Cube.zcoord[indexz[0]]
        xloc = Cube.xcoord[indexx[0]]
        yloc = Cube.ycoord[indexy[0]]


        v2ab_transform = Cube.transform[ifile]

        if(self.coord_system == 'alpha-beta'):
            distance1 = abs(xloc - alpha)
            distance2 = abs(yloc - beta)

        if(self.coord_system == 'v2-v3'):
            distance1 = abs(xloc - coord1)
            distance2 = abs(yloc - coord2)

        distance3 = abs(zloc - wave)
        distance12 = distance1 * distance1
        distance22 = distance2 * distance2
        distance32 = distance3 * distance3

        iz = 0
        for zz in indexz[0]:
            istart = zz * nplane
            iy = 0
            for yy in indexy[0]:
                ix = 0
                for xx in indexx[0]:

#                    if(self.coord_system=='v2-v3'):
#                        xloc60 = xloc[ix]/60.0
#                        yloc60 = yloc[iy]/60.0
#                        alpha_spaxel,beta_spaxel,wave_spaxel = v2ab_transform(xloc60,yloc60,zloc[iz])
#                        alpha_distance = abs(alpha-alpha_spaxel)
#                        beta_distance = abs(beta-beta_spaxel)
#                        wave_distance  = abs(wave-wave_spaxel)

#                        xn = alpha_distance/weight_alpha
#                        yn = beta_distance/weight_beta
#                        wn = wave_distance/weight_wave
#                        weight_distance = xn*xn + yn*yn + wn*wn


                    weight_distance = distance12[ix] + distance22[iy] + distance32[iz]

                    if(weight_distance < lower_limit): weight_distance = lower_limit
                    weight_distance = 1.0 / weight_distance


                    cube_index = istart + yy * Cube.naxis1 + xx
                    spaxel[cube_index].ipointcloud.append(ipt)
                    spaxel[cube_index].pointcloud_weight.append(weight_distance)

                    #if(zz == 100 and xx == 7 and yy == 10):
                    #    print('Match',ipt,cube_index)
                    #    print('Pt point',coord1,coord2,wave)
                    #    print('Pt a-b-l',alpha,beta)
                    #    print('x y',x+1,y+1)
                    #    print(xloc[ix],xloc[ix]-coord1)
                    #    print(yloc[iy],yloc[iy]-coord2)
                    #    print(zloc[iz],zloc[iz]-wave)

                    ix = ix + 1
                    iprint = iprint + 1
#                    if(iprint > 10000):
#                        iprint = 0
                iy = iy + 1
            iz = iz + 1

#_______________________________________________________________________

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

