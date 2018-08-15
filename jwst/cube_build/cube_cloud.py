
 # Routines used in Spectral Cube Building

#import sys
import numpy as np
import logging
from ..assign_wcs import nirspec
from ..datamodels import dqflags
from . import coord
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
#________________________________________________________________________________
#********************************************************************************
# HELPER ROUTINES for CubeData class defined in cube_build.py
# these methods relate to wcs type procedures.
# determine_scale
#********************************************************************************
def match_det2cube(self, input_model,
                  x, y, file_slice_no,
                  this_par1, this_par2,
                  spaxel,
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
    if self.instrument == 'MIRI':
        coord1, coord2, wave = input_model.meta.wcs(x, y) # for entire detector find  ra,dec,lambda
        valid1 = ~np.isnan(coord1)

# MIRI PSF weighting
# weighting based on distance in alpha-beta plane
        if self.weighting == 'miripsf':
            # get transforms needed for this weighting

            det2ab_transform = input_model.meta.wcs.get_transform('detector', 'alpha_beta')
            worldtov23 = input_model.meta.wcs.get_transform("world", "v2v3")
            v2ab_transform = input_model.meta.wcs.get_transform('v2v3',
                                                                'alpha_beta')
            coord1, coord2, wave = det2ab_transform(x, y)
            valid1 = ~np.isnan(coord1)
            wave_resol = self.instrument_info.Get_RP_ave_Wave(this_par1, this_par2)
            alpha_resol = self.instrument_info.Get_psf_alpha_parameters()
            beta_resol = self.instrument_info.Get_psf_beta_parameters()

            # transform Cube Spaxel centers to alpha,beta of exposure
            # for MIRI weighting parameters are based on distance in
            # alpha-beta coord system
            # transform the cube coordinate values to alpha and beta values
            # xi,eta -> ra,dec
            # ra-dec -> v2,v3
            # v2,v3 -> local alpha,beta

        if self.coord_system == 'alpha-beta': # making Alpha-Beta final Cubes
            det2ab_transform = input_model.meta.wcs.get_transform('detector',
                                                                  'alpha_beta')
            coord1, coord2, wave = det2ab_transform(x, y)
            valid1 = ~np.isnan(coord1)
#________________________________________________________________________________
    elif self.instrument == 'NIRSPEC': # for NIRSPEC need to look over every slice
        # Slice No are provided by ifu_cube.py
        islice = file_slice_no
        slice_wcs = nirspec.nrs_wcs_set_input(input_model, islice)

        coord1, coord2, wave = slice_wcs(x, y) # return v2,v3 are in degrees
        valid1 = ~np.isnan(coord1)

#________________________________________________________________________________
#________________________________________________________________________________
# Slices are curved on detector. A slice region is grabbed by corner regions so
# the region returned may include pixels not value for slice. There are gaps
# between the slices. Pixels not belonging to a slice are assigned NaN values.

    x = x[valid1]
    y = y[valid1]
    coord1 = coord1[valid1]
    coord2 = coord2[valid1]
    wave = wave[valid1]

#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded
# from the cube mapping

    flux_all = input_model.data[y, x]
    dq_all = input_model.dq[y, x]
    valid4 = np.isfinite(flux_all)

    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] +
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] +
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])

    # good data holds the location of pixels we want to map to cube
    good_data = np.where((np.bitwise_and(dq_all, all_flags) == 0) & (valid4))
#________________________________________________________________________________
# Now get the final values to use for cube_building
    x = x[good_data]
    y = y[good_data]
    coord1 = coord1[good_data]
    coord2 = coord2[good_data]
    wave = wave[good_data]

    flux = input_model.data[y, x]
#    error = input_model.err[y, x]  - we might use this later

# for now comment out - check if we are ever going to use c1_offset, c2_offset
#    ra = ra - c1_offset/3600.0
#    dec = dec - c2_offset/3600.0

    # find the Cube Coordinates
    if self.coord_system == 'ra-dec':
# xi,eta in arc seconds (works both for MIRI and NIRSPEC)
        # need to redefine what coord1 and coord2 is if working
        # with Ra,Dec (majority of cases)
        xi, eta = coord.radec2std(self.Crval1, self.Crval2, coord1, coord2)
        coord1 = xi
        coord2 = eta

    nplane = self.naxis1 * self.naxis2
    lower_limit = 0.01

#    iprint = 0
# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn = coord1.size
#    print('looping over n points on exposure and mapping to cloud',nn)
#________________________________________________________________________________
    for ipt in range(0, nn - 1):
#________________________________________________________________________________
        # Cube.Xcenters, Ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.
        # Find the spaxels that fall withing ROI of point cloud defined  by
        # coord1,coord2,wave

        xdistance = (self.Xcenters - coord1[ipt])
        ydistance = (self.Ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        #indexr holds the index of the Xcenters,Ycenters points that fall in ROI
        #indexz holds the index of the cube wavelength (zcoord) points that
        # fall in the ROI
        indexr = np.where(radius <= self.rois)
        indexz = np.where(abs(self.zcoord - wave[ipt]) <= self.roiw)

        # Pull out the Spaxel center values that are used for weighting
        zlam = self.zcoord[indexz]
        xi_cube = self.Xcenters[indexr]
        eta_cube = self.Ycenters[indexr]

        # Form the distances that are used for weighting the spaxel points
        d1 = (xi_cube - coord1[ipt]) / self.Cdelt1
        d2 = (eta_cube - coord2[ipt]) / self.Cdelt2
        d3 = (zlam - wave[ipt]) / self.Cdelt3

        dxy = d1 * d1 + d2 * d2

        # Form arrays to calculate weight factor
        dxy_matrix = np.tile(dxy[np.newaxis].T, [1, d3.shape[0]])
        d3_matrix = np.tile(d3 * d3, [dxy_matrix.shape[0], 1])

        wdistance = dxy_matrix + d3_matrix
        weight_distance = np.power(np.sqrt(wdistance), self.weight_power)
        weight_distance[weight_distance < lower_limit] = lower_limit
        weight_distance = 1.0 / weight_distance

        izindex = (indexz[0] * nplane).astype(np.int)
        idxy_matrix = np.tile(indexr[0][np.newaxis].T, [1, d3.shape[0]])
        id3_matrix = np.tile(izindex, [dxy_matrix.shape[0], 1])
        icube_index = idxy_matrix + id3_matrix

        weight_distance = weight_distance.flatten()
        icube_index = (icube_index.flatten()).astype(np.int)
#**********************************************************************
# WEIGHTING = MIRPSF
#TODO with the CPD 7 delivery of new reference file is this going away ???
# if it stays then update to use list comprehension or array math
#**********************************************************************
        if self.weighting == 'miripsf':
            weights = FindNormalizationWeights(wave[ipt],
                                               wave_resol,
                                               alpha_resol,
                                               beta_resol)
            weight_alpha, weight_beta, weight_wave = weights
#________________________________________________________________________________
# loop over the points in the ROI
            for iz, zz in enumerate(indexz[0]):
                weight_distance = []
#                istart = zz * nplane
                for ir, rr in enumerate(indexr[0]):
#________________________________________________________________________________
# if weight is miripsf -distances determined in alpha-beta coordinate system
                    ra_spaxel, dec_spaxel = coord.std2radec(self.Crval1,
                                                            self.Crval2,
                                                            xi_cube[ir],
                                                            eta_cube[ir])

                    v2_spaxel, v3_spaxel, zl = worldtov23(ra_spaxel,
                                                        dec_spaxel,
                                                        zlam[iz])

                    alpha_spaxel, beta_spaxel, wave_spaxel = v2ab_transform(v2_spaxel,
                                                                        v3_spaxel,
                                                                        zlam[iz])
                    alpha_distance = coord1[ipt] - alpha_spaxel
                    beta_distance = coord2[ipt] - beta_spaxel
                    wave_distance = abs(wave[ipt] - wave_spaxel)

                    xn = alpha_distance / weight_alpha
                    yn = beta_distance / weight_beta
                    wn = wave_distance / weight_wave

                    this_wdistance = xn * xn + yn * yn + wn * wn
                    weight_distance.append(this_wdistance)
#________________________________________________________________________________
# If debug information is wanted find cube Location for later
# comment out now for speed
#________________________________________________________________________________
#        if self.debug_pixel == 1:
#            yy_cube = (indexr[0]/self.naxis1).astype(np.int)
#            xx_cube = indexr[0] -yy_cube*self.naxis1
#            cube_location = [ (xx,yy,zz) for zz in indexz[0] for xx,yy in zip(xx_cube,yy_cube)]

#----------------------------------------------------------------------
        # loop over all spaxels that overlap with ipt
        # TODO: can this be done faster using some other method than looping
        
        for iz, zz in enumerate(icube_index):
            spaxel[zz].flux = spaxel[zz].flux + weight_distance[iz] * flux[ipt]
            spaxel[zz].flux_weight = spaxel[zz].flux_weight + weight_distance[iz]
            spaxel[zz].iflux = spaxel[zz].iflux + 1
#----------------------------------------------------------------------
# debug informatons - comment out now for speed
#            if( self.debug_pixel == 1 and self.xdebug == cube_location[iz][0] and
#                self.ydebug == cube_location[iz][1] and
#                self.zdebug ==cube_location[iz][2]):

#                log.debug('For spaxel %d %d %d, %d detector x,y,flux %d %d %f %d %f '
#                          %(self.xdebug+1,self.ydebug+1,
#                            self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                            flux[ipt],file_slice_no,weight_distance[iz]))
#                self.spaxel_debug.write(
#                    'For spaxel %d %d %d %d, detector x,y,flux %d %d %f %d %f %f %f %f %f '
#                    %(self.xdebug+1,self.ydebug+1,
#                      self.zdebug+1,ipt,xpix[ipt]+1,ypix[ipt]+1,
#                      flux[ipt],file_slice_no,weight_distance[iz],wave[ipt],
#                      d1*self.Cdelt1,d2*self.Cdelt2,d3*self.Cdelt3) +' \n')

#----------------------------------------------------------------------
# comment out for speed
#        iprint = iprint +1
#        if iprint == 10000:
#            log.debug('Mapping point and finding ROI for point cloud # %d %d' %(ipt,nn))
#            iprint = 0
#_______________________________________________________________________
def FindWaveWeights(channel, subchannel):
    """
    Short Summary
    -------------
    Get the wavelength normalization weights that we will
    use to normalize wavelengths.

    Parameters
    ----------
    channel- channel for point
    subchannel- subchannel for point

    Returns
    -------
    normalized weighting for wavelength for this channel, subchannel

    """

    if channel == '1':
        if subchannel == 'SHORT':
            a = 3050.0
            c = 3340.0
            wa = 4.91
            wc = 5.79
        elif subchannel == 'MEDIUM':
            a = 2920.0
            c = 3400.0
            wa = 5.6
            wc = 6.62
        elif subchannel == 'LONG':
            a = 2800.0
            c = 3220.0
            wa = 6.46
            wc = 7.63


    if channel == '2':
        if subchannel == 'SHORT':
            a = 2700.0
            c = 2800.0
            wa = 7.55
            wc = 8.91
        elif subchannel == 'MEDIUM':
            a = 2600.0
            c = 2880.0
            wa = 8.71
            wc = 10.34
        elif subchannel == 'LONG':
            a = 2590.0
            c = 3000.0
            wa = 9.89
            wc = 11.71

    if channel == '3':
        if subchannel == 'SHORT':
            a = 2390.0
            c = 2650.0
            wa = 11.50
            wc = 13.59
        elif subchannel == 'MEDIUM':
            a = 1600.0
            c = 2400.0
            wa = 13.19
            wc = 15.58
        elif subchannel == 'LONG':
            a = 1850.0
            c = 2550.0
            wa = 15.40
            wc = 18.14

    if channel == '4':
        if subchannel == 'SHORT':
            a = 1320.0
            c = 1720.0
            wa = 17.88
            wc = 21.34
        elif subchannel == 'MEDIUM':
            a = 1550.0
            c = 1600.0
            wa = 20.69
            wc = 24.68
        elif subchannel == 'LONG':
            a = 1450.0
            c = 1200.0
            wa = 23.83
            wc = 28.43
    return a, c, wa, wc

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

    # alpha psf weighting
    alpha_wave_cutoff = alpha_resol[0]
    alpha_a_short = alpha_resol[1]
    alpha_b_short = alpha_resol[2]
    alpha_a_long = alpha_resol[3]
    alpha_b_long = alpha_resol[4]
    if wavelength < alpha_wave_cutoff:
        alpha_weight = alpha_a_short + (alpha_b_short * wavelength)
    else:
        alpha_weight = alpha_a_long + (alpha_b_long * wavelength)

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
    resolution = a_ave + (b_ave * wave_diff) + (c_ave * wave_diff * wave_diff)
    lambda_weight = wavelength / resolution

    weight = [alpha_weight, beta_weight, lambda_weight]
    return weight
