
 # Routines used in Spectral Cube Building
import numpy as np
import math
import logging
from .. import datamodels
from ..assign_wcs import nirspec
from ..datamodels import dqflags
from . import coord
from gwcs import wcstools

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
    file_slice_no: the index on the files that are used to construct the Cube
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

        det2ab_transform = input_model.meta.wcs.get_transform('detector', 'alpha_beta')
        detector2v23 = input_model.meta.wcs.get_transform('detector', 'v2v3')
        v23toworld = input_model.meta.wcs.get_transform("v2v3", "world")
        worldtov23 = input_model.meta.wcs.get_transform("world", "v2v3")
        v2ab_transform = input_model.meta.wcs.get_transform('v2v3',
                                                            'alpha_beta')

        alpha, beta, wave = det2ab_transform(x, y)
        v2, v3, lam23 = detector2v23(x, y)
        ra, dec, lam = v23toworld(v2, v3, lam23)

        valid1 = np.isfinite(v2)
        valid2 = np.isfinite(v3)

        if self.weighting == 'miripsf':
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

    elif self.instrument == 'NIRSPEC':
        islice = file_slice_no
        slice_wcs = nirspec.nrs_wcs_set_input(input_model, islice)

#        x, y = wcstools.grid_from_bounding_box(slice_wcs.bounding_box,step=(1,1),center=True)


        ra, dec, lam = slice_wcs(x, y) # return v2,v3 are in degrees
        valid1 = np.isfinite(ra)
        valid2 = np.isfinite(dec)
#________________________________________________________________________________
#________________________________________________________________________________
# Slices are curved on detector. A slice region is grabbed by corner regions so
# the region returned may include pixels not value for slice. There are gaps
# between the slices. Pixels not belonging to a slice are assigned NaN values.

    x = x.astype(np.int)
    y = y.astype(np.int)

    flux_all = input_model.data[y, x]
#    error_all = input_model.err[y, x]
    dq_all = input_model.dq[y, x]

    valid3 = np.isfinite(lam)
    valid4 = np.isfinite(flux_all)
    valid = valid1 & valid2 & valid3 & valid4
#________________________________________________________________________________
# using the DQFlags from the input_image find pixels that should be excluded
# from the cube mapping
    all_flags = (dqflags.pixel['DO_NOT_USE'] + dqflags.pixel['DROPOUT'] +
                 dqflags.pixel['NON_SCIENCE'] +
                 dqflags.pixel['DEAD'] + dqflags.pixel['HOT'] +
                 dqflags.pixel['RC'] + dqflags.pixel['NONLINEAR'])

    # find the location of all the values to reject in cube building
    good_data = np.where((np.bitwise_and(dq_all, all_flags) == 0) & (valid == True))

    # good data holds the location of pixels we want to map to cube
    flux = flux_all[good_data]
#    error = error_all[good_data]
    wave = lam[good_data]

    xpix = x[good_data] # only used for testing
    ypix = y[good_data] # only used for testing

    ra = ra - c1_offset / 3600.0
    dec = dec - c2_offset / 3600.0
    ra_use = ra[good_data]
    dec_use = dec[good_data]
    if self.instrument == 'MIRI':
        # need alpha,beta if weigthing is miripsf or cubes in alpha-beta space
        alpha_det = alpha[good_data]
        beta_det = beta[good_data]
# MIRI can make cubes in alpha-beta:
    if self.coord_system == 'alpha-beta':
        coord1 = alpha[good_data]
        coord2 = beta[good_data]

    else:
# xi,eta in arc seconds
        xi, eta = coord.radec2std(self.Crval1, self.Crval2, ra_use, dec_use)
        coord1 = xi
        coord2 = eta

    nplane = self.naxis1 * self.naxis2
    lower_limit = 0.01

    iprint = 0
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
        xdistance = (self.Xcenters - coord1[ipt])
        ydistance = (self.Ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        indexr = np.where(radius <= self.rois)
        indexz = np.where(abs(self.zcoord - wave[ipt]) <= self.roiw)

#        print('indexz',indexz)
#        print('indexr',indexr)
        zlam = self.zcoord[indexz]        # z Cube values falling in wavelength roi
        xi_cube = self.Xcenters[indexr]   # x Cube values within radius
        eta_cube = self.Ycenters[indexr]  # y cube values with the radius

#        print('found xi_cube',xi_cube)
#        print('found eta_cube',eta_cube)

#________________________________________________________________________________
# loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for ir, rr in enumerate(indexr[0]):
                yy_cube = int(rr / self.naxis1)
                xx_cube = rr - yy_cube * self.naxis1
#                print('xx yy cube',rr,self.naxis1,xx_cube,yy_cube)
#________________________________________________________________________________
                if self.weighting == 'msm':
                    d1 = (xi_cube[ir] - coord1[ipt]) / self.Cdelt1
                    d2 = (eta_cube[ir] - coord2[ipt]) / self.Cdelt2
                    d3 = (zlam[iz] - wave[ipt]) / self.Cdelt3

                    weight_distance = math.sqrt(d1 * d1 + d2 * d2 + d3 * d3)
                    weight_distance = math.pow(weight_distance, self.weight_power)
#________________________________________________________________________________
# if weight is miripsf -distances determined in alpha-beta coordinate system
                elif self.weighting == 'miripsf':
                    weights = FindNormalizationWeights(wave[ipt],
                                                       wave_resol,
                                                       alpha_resol,
                                                       beta_resol)


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
                    alpha_distance = alpha_det[ipt] - alpha_spaxel
                    beta_distance = beta_det[ipt] - beta_spaxel
                    wave_distance = abs(wave[ipt] - wave_spaxel)

                    xn = alpha_distance / weights[0]
                    yn = beta_distance / weights[1]
                    wn = wave_distance / weights[2]

                        # only included the spatial dimensions
                    weight_distance = math.sqrt(xn * xn + yn * yn + wn * wn)
                    weight_distance = math.pow(weight_distance, self.weight_power)
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
