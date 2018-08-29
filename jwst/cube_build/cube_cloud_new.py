 # Routines used in Spectral Cube Building
import numpy as np
import math
import logging
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
#________________________________________________________________________________
def match_det2cube_msm(naxis1,naxis2,naxis3,
                       cdelt1,cdelt2,cdelt3,
                       rois,roiw,
                       msm_weight_power,
                       xcenters,ycenters,zcoord,
                       spaxel,flux,
                       coord1,coord2,wave):
                       

    """
    Short Summary
    -------------
    Match the Point Cloud members to the spaxel centers that fall in the ROI.
    For each spaxel the coord1,coord1 and wave point cloud members are weighting
    according to modified shepard method of inverse weighting based on the distance between
    the point cloud member and the spaxel center. 

    Parameters
    ----------
    naxis1,naxis2,naxis3: size of the ifucube
    cdelt1,cdelt2,cdelt3: ifucube spaxel size in the 3 dimensions
    rois, roiw: region of influence size in spatial and spectral dimension
    weight_power: msm weighting parameter
    xcenter,ycenter: spaxel center locations  in 1st and 2nd dimensions. These values are 2 X 2 grid
    spaxel centers. 
    zcoord: spaxel center locations in 3rd dimensions
    spaxel: class holding ifucube spaxel values

    Returns
    -------
    spaxel class matched to detector pixels with flux and weighting updated for each
    match


    """

    nplane = naxis1 * naxis2
    lower_limit = 0.01

# loop over the spatial dimensions of the ifucube first

#    iprint = 0
# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn = coord1.size

    print('looping over n points',nn)
#________________________________________________________________________________
    for ipt in range(0, nn - 1):
#________________________________________________________________________________
        # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.
        # find the point cloud members falling on within ROIS of this spaxel
        # The coordinates of the point cloud are defined by coord1,coord2,wave

        xdistance = (xcenters[ispaxel] - coord1)
        ydistance = (ycenters[ispaxel] - coord2)
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        indexr = np.where(radius <= rois)

        # 
        indexz = np.where(abs(zcoord[ispaxel] - wave) <= roiw)

#        print('indexz',indexz)
#        print('indexr',indexr)
        wave_found = wave[indexz]        # z Cube values falling in wavelength roi
        xi_found = coord1[indexr]   # x Cube values within radius
        eta_found = coord2[indexr]  # y cube values with the radius
#________________________________________________________________________________
        # form the arrays to be used in calculating the weighting. 
        
        d1 = np.array(xi_found - xcenter[ispaxel]) / cdelt1
        d2 = np.array(eta_found - ycenter[ispaxel]) / cdelt2
        d3 = np.array(wave_found - wave[ispaxel]) / cdelt3

        dxy = (d1 * d1) + (d2 * d2)
        dxy_matrix = np.tile(dxy[np.newaxis].T,[1,d3.shape[0]])
        d3_matrix = np.tile(d3 * d3,[dxy_matrix.shape[0],1])

        wdistance = dxy_matrix + d3_matrix
        weight_distance = np.power(np.sqrt(wdistance), msm_weight_power)
        weight_distance[weight_distance < lower_limit] = lower_limit
        weight_distance = 1.0 / weight_distance


        # determine the spaxel xx_cube,yy_cube values of these spaxels in
        # the ROI so they can be used to pull out the flux of the median
        # sky cube.
                yy_cube = (indexr[0] / self.naxis1).astype(np.int)
                xx_cube = indexr[0] - yy_cube * self.naxis1
                scf = np.array([self.cube_flux[zz, yy_cube[ir], xx_cube[ir]]
                                for ir, rr in enumerate(indexr[0]) for zz in indexz[0]])
                scf = np.reshape(scf, weight_distance.shape)

#________________________________________________________________________________
# loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for ir, rr in enumerate(indexr[0]):
#                yy_cube = int(rr / self.naxis1)
#                xx_cube = rr - yy_cube * self.naxis1
#                print('xx yy cube',rr,self.naxis1,xx_cube,yy_cube)
#________________________________________________________________________________

#________________________________________________________________________________
# We have found the weight_distance based on instrument type
                if weight_distance < lower_limit: weight_distance = lower_limit
                weight_distance = 1.0 / weight_distance

                cube_index = istart + rr
                spaxel[cube_index].flux = spaxel[cube_index].flux + weight_distance * flux[ipt]
                spaxel[cube_index].flux_weight = spaxel[cube_index].flux_weight + weight_distance
                spaxel[cube_index].iflux = spaxel[cube_index].iflux + 1
#                if( self.debug_pixel == 1 and self.xdebug == xx_cube and
#                      self.ydebug == yy_cube and self.zdebug == zz):

#                    log.debug('For spaxel %d %d %d, %d detector x,y,flux %d %d %f %d %f '
#                                  %(self.xdebug+1,self.ydebug+1,
#                                    self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                                    flux[ipt],file_slice_no,weight_distance))
#                    self.spaxel_debug.write('For spaxel %d %d %d %d, detector x,y,flux %d %d %f %d %f %f %f %f %f %f '
#                                  %(self.xdebug+1,self.ydebug+1,
#                                    self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                                    flux[ipt],file_slice_no,weight_distance,wave[ipt],zlam[iz],
#                                    d1*self.Cdelt1,d2*self.Cdelt2,d3*self.Cdelt3) +' \n')
#        iprint = iprint + 1
#        if iprint == 10000:
#            log.debug('Finding ROI for point cloud # %d %d' % (ipt , nn))
#            iprint = 0
#_______________________________________________________________________

def match_det2cube_miripsf(alpha_resol,beta_resol,wave_resol,
                           worldtov23,
                           v2ab_tranform,
                           naxis1,naxis2,naxis3,
                           cdelt1,cdelt2,cdelt3,
                           crval1,crval2,
                           rois,roiw,
                           weight_power,
                           xcenters,ycenters,zcoord,
                           spaxel,
                           coord1,coord2,wave,alpha_det,beta_det):
    """
    Short Summary

    -------------
    Map coordinates coord1,coord2, and wave of the point cloud to which spaxels they
    overlap with in the ifucube.
    For each spaxel the coord1,coord1 and wave point cloud members are weighting
    according to the miri psf and lsf. 
    The weighting function is based on the distance the point cloud member and spaxel
    center in the alph-beta coordinate system.
    The alpha and beta value of each point cloud member is passed to this routine
    The alpha and beta value of the spaxel center is determined from the passed
    in transforms: worldtov23 and v2ab_transform 

    Parameters
    ----------
    alpha_resol,beta_resol,wave_resol: alpha,beta and wavelength resolution table
    worldtov23: transform ra, dec -> v2,v3
    v2ab_tranform: transform v2,v3 -> alpha,beta on miri detector plane
    naxis1,naxis2,naxis3: size of the ifucube
    cdelt1,cdelt2,cdelt3: ifucube spaxel size in the 3 dimensions
    crval1,crval2: ra and dec center of ifu cube used to transform xi,eta spaxel -> ra,dec
    rois, roiw: region of influence size in spatial and spectral dimension
    weight_power: msm weighting parameter
    xcenter,ycenter: spaxel center locations in 1st and 2nd dimension. These values are 2 X2 grid
    of spaxel center locations.
    zcoord: spaxel center locations in 3rd dimension
    spaxel: class holding ifucube spaxel values
    coord1,coord2,wave pixel coordinates mapped to output frame
    alpha_det,beta_det alpha,beta values of pixel coordinates:

    Returns
    -------
    spaxel class matched to detector pixels with flux and weighting updated for each
    match

    """

    nplane = naxis1 * naxis2
    lower_limit = 0.01

#    iprint = 0
# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn = coord1.size

#    print('looping over n points mapping to cloud',nn)
#________________________________________________________________________________
    for ipt in range(0, nn - 1):
#________________________________________________________________________________
        # Cube.Xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.
        # find the spaxels that fall withing ROI of point cloud defined  by
        # coord1,coord2,wave
#        if(ipt > 2): sys.exit('STOP')
#        print('For point ',coord1[ipt],coord2[ipt],wave[ipt],ipt)

#        if(ipt == 0):
#            print('size of Xcenters',self.Xcenters.size)
        xdistance = (xcenters - coord1[ipt])
        ydistance = (ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        indexr = np.where(radius <= rois)
        indexz = np.where(abs(zcoord - wave[ipt]) <= roiw)

#        print('indexz',indexz)
#        print('indexr',indexr)
        zlam = zcoord[indexz]        # z Cube values falling in wavelength roi
        xi_cube = xcenters[indexr]   # x Cube values within radius
        eta_cube = ycenters[indexr]  # y cube values with the radius

#        print('found xi_cube',xi_cube)
#        print('found eta_cube',eta_cube)

#________________________________________________________________________________
# loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for ir, rr in enumerate(indexr[0]):
#                yy_cube = int(rr / self.naxis1)
#                xx_cube = rr - yy_cube * self.naxis1
#                print('xx yy cube',rr,self.naxis1,xx_cube,yy_cube)
#________________________________________________________________________________
#________________________________________________________________________________
# if weight is miripsf -distances determined in alpha-beta coordinate system

                weights = FindNormalizationWeights(wave[ipt],
                                                   wave_resol,
                                                   alpha_resol,
                                                   beta_resol)

                
                ra_spaxel, dec_spaxel = coord.std2radec(crval1,
                                                        crval2,
                                                        xi_cube[ir],
                                                        eta_cube[ir])

                v2_spaxel, v3_spaxel, zl = worldtov23(ra_spaxel,
                                                      dec_spaxel,
                                                      zlam[iz])
                
                alpha_spaxel, beta_spaxel, wave_spaxel = v2ab_transform(v2_spaxel,
                                                                        v3_spaxel,
                                                                        zlam[iz])
                alpha_distance = alpha_det[ipt] - alpha_spaxel
                beta_distance = beta_det[ipt] - beta_spaxel
                wave_distance = abs(wave[ipt] - wave_spaxel)

                xn = alpha_distance / weights[0]
                yn = beta_distance / weights[1]
                wn = wave_distance / weights[2]

                        # only included the spatial dimensions
                weight_distance = math.sqrt(xn * xn + yn * yn + wn * wn)
                weight_distance = math.pow(weight_distance, weight_power)
#________________________________________________________________________________
# We have found the weight_distance based on instrument type

                if weight_distance < lower_limit: weight_distance = lower_limit
                weight_distance = 1.0 / weight_distance

                cube_index = istart + rr
                spaxel[cube_index].flux = spaxel[cube_index].flux + weight_distance * flux[ipt]
                spaxel[cube_index].flux_weight = spaxel[cube_index].flux_weight + weight_distance
                spaxel[cube_index].iflux = spaxel[cube_index].iflux + 1


#                if( self.debug_pixel == 1 and self.xdebug == xx_cube and
#                      self.ydebug == yy_cube and self.zdebug == zz):

#                    log.debug('For spaxel %d %d %d, %d detector x,y,flux %d %d %f %d %f '
#                                  %(self.xdebug+1,self.ydebug+1,
#                                    self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                                    flux[ipt],file_slice_no,weight_distance))
#                    self.spaxel_debug.write('For spaxel %d %d %d %d, detector x,y,flux %d %d %f %d %f %f %f %f %f %f '
#                                  %(self.xdebug+1,self.ydebug+1,
#                                    self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                                    flux[ipt],file_slice_no,weight_distance,wave[ipt],zlam[iz],
#                                    d1*self.Cdelt1,d2*self.Cdelt2,d3*self.Cdelt3) +' \n')
#        iprint = iprint + 1
#        if iprint == 10000:
#            log.debug('Finding ROI for point cloud # %d %d' % (ipt , nn))
#            iprint = 0
#_______________________________________________________________________
def FindNormalizationWeights(wavelength,
                              wave_resol,
                              alpha_resol,
                              beta_resol):
    """
    Short Summary
    -------------
    we need to normalize how to each point cloud in the spaxel. The
    normalization of weighting is determined from width of PSF as well
    as wavelength resolution

    Parameters
    ----------
    wave_resol:
    alpha_resol:
    beta_resol:

    Returns
    -------
    normalized weighting for 3 dimension

    """
    alpha_weight = 1.0
    beta_weight = 1.0
    lambda_weight = 1.0

    # alpha psf weighting
    alpha_wave_cutoff = alpha_resol[0]
    alpha_a_short = alpha_resol[1]
    alpha_b_short = alpha_resol[2]
    alpha_a_long = alpha_resol[3]
    alpha_b_long = alpha_resol[4]
    if wavelength < alpha_wave_cutoff:
        alpha_weight = alpha_a_short + alpha_b_short * wavelength
    else:
        alpha_weight = alpha_a_long + alpha_b_long * wavelength

    # beta psf weighting
    beta_wave_cutoff = beta_resol[0]
    beta_a_short = beta_resol[1]
    beta_b_short = beta_resol[2]
    beta_a_long = beta_resol[3]
    beta_b_long = beta_resol[4]
    if wavelength < beta_wave_cutoff:
        beta_weight = beta_a_short + beta_b_short * wavelength
    else:
        beta_weight = beta_a_long + beta_b_long * wavelength

    # wavelength weighting
    wavecenter = wave_resol[0]
    a_ave = wave_resol[1]
    b_ave = wave_resol[2]
    c_ave = wave_resol[3]
    wave_diff = wavelength - wavecenter
    resolution = a_ave + b_ave * wave_diff + c_ave * wave_diff * wave_diff
    lambda_weight = wavelength / resolution
    weight = [alpha_weight, beta_weight, lambda_weight]
    return weight
