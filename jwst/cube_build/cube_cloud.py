""" Map the detector pixels to the cube coordinate system
"""
import numpy as np
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def match_det2cube_msm(naxis1, naxis2, naxis3,
                       cdelt1, cdelt2,
                       zcdelt3,
                       xcenters, ycenters, zcoord,
                       spaxel_flux,
                       spaxel_weight,
                       spaxel_iflux,
                       flux,
                       coord1, coord2, wave,
                       rois_pixel, roiw_pixel, weight_pixel, softrad_pixel):

    """ Map the detector pixels to the cube spaxels using the MSM parameters


    Match the Point Cloud members to the spaxel centers that fall in the ROI.
    For each spaxel the coord1,coord1 and wave point cloud members are weighed
    according to modified shepard method of inverse weighting based on the
    distance between the point cloud member and the spaxel center.

    Parameters
    ----------
    naxis1 : int
       size of the ifucube in 1st axis
    naxis2 : int
       size of the ifucube in 2nd axis
    naxis3 : int
       size of the ifucube in 3rd axis
    cdelt1 : float
       ifucube spaxel size in axis 1 dimension
    cdelt2 : float
       ifucube spaxel size in axis 2 dimension
    cdelt3_normal :float
       ifu spectral size at wavelength
    rois_pixel : float
       region of influence size in spatial dimension
    roiw_pixel : float
       region of influence size in spectral dimension
    weight_power : float
       msm weighting parameter
    xcenter : numpy.ndarray
       spaxel center locations 1st dimensions.
    ycenter : numpy.ndarray
       spaxel center locations 2nd dimensions.
    zcoord : numpy.ndarray
        spaxel center locations in 3rd dimensions
    spaxel_flux : numpy.ndarray
       contains the weighted summed detector fluxes that fall
       within the roi
    spaxel_weight : numpy.ndarray
       contains the summed weights assocated with the detector fluxes
    spaxel_iflux : numpy.ndarray
       number of detector pixels falling with roi of spaxel center
    flux : numpy.ndarray
       array of detector fluxes associated with each position in
       coorr1, coord2, wave
    coord1 : numpy.ndarray
       contains the spatial coordinate for 1st dimension for the mapped
       detector pixel
    coord2 : numpy.ndarray
       contains the spatial coordinate for 2nd dimension for the mapped
       detector pixel
    wave : numpy.ndarray
       contains the spectral coordinate  for the mapped detector pixel

    Returns
    -------
    spaxel_flux, spaxel_weight, and spaxel_ifux updated with the information
    from the detector pixels that fall within the roi if the spaxel center.
    """

    nplane = naxis1 * naxis2

    # now loop over the pixel values for this region and find the spaxels that fall
    # within the region of interest.
    nn = coord1.size
    # Left this commented code in to check when we get new cube par reference files
    # The roi size was too small at long wavelength that may pixels did not match
    # to a spaxel.
    #    ilow = 0
    #    ihigh = 0
    #    imatch  = 0
    #    print('looping over n points mapping to cloud',nn)
    # ________________________________________________________________________________
    for ipt in range(0, nn - 1):
        # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.
        # find the spaxels that fall withing ROI of point cloud defined  by
        # coord1,coord2,wave
        lower_limit = softrad_pixel[ipt]
        xdistance = (xcenters - coord1[ipt])
        ydistance = (ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        indexr = np.where(radius <= rois_pixel[ipt])
        indexz = np.where(abs(zcoord - wave[ipt]) <= roiw_pixel[ipt])

        # on the wavelength boundaries the point cloud may not be in the IFUCube
        # the edge cases are skipped and not included in final IFUcube to avoid noisy results
        # Left commented code for checking later for NIRSPEC the spectral size
        # in the reference file may be too small
#        if len(indexz[0]) == 0:
#            if wave[ipt] < zcoord[0]:
#               ilow = ilow + 1
#            elif wave[ipt] > zcoord[-1]:
#               ihigh = ihigh + 1
#            else:
#                imatch = imatch + 1
#                print(' no z match found ',wave[ipt],roiw_pixel[ipt])
#                print(zcoord[naxis3-11:naxis3])
#                diff = abs(zcoord[naxis3-11:naxis3] - wave[ipt])
#                print(diff)
#                exit()
        if len(indexz[0]) > 0:
            d1 = np.array(coord1[ipt] - xcenters[indexr]) / cdelt1
            d2 = np.array(coord2[ipt] - ycenters[indexr]) / cdelt2
            d3 = np.array(wave[ipt] - zcoord[indexz]) / zcdelt3[indexz]

            dxy = (d1 * d1) + (d2 * d2)

        # shape of dxy is #indexr or number of overlaps in spatial plane
        # shape of d3 is #indexz or number of overlaps in spectral plane
        # shape of dxy_matrix & d3_matrix  (#indexr, #indexz)
        # rows = number of overlaps in spatial plane
        # cols = number of overlaps in spectral plane
            dxy_matrix = np.tile(dxy[np.newaxis].T, [1, d3.shape[0]])
            d3_matrix = np.tile(d3 * d3, [dxy_matrix.shape[0], 1])

            wdistance = dxy_matrix + d3_matrix
            weight_distance = np.power(np.sqrt(wdistance), weight_pixel[ipt])
            weight_distance[weight_distance < lower_limit] = lower_limit
            weight_distance = 1.0 / weight_distance
            weight_distance = weight_distance.flatten('F')
            weighted_flux = weight_distance * flux[ipt]

            icube_index = [iz * nplane + ir for iz in indexz[0] for ir in indexr[0]]
            spaxel_flux[icube_index] = spaxel_flux[icube_index] + weighted_flux
            spaxel_weight[icube_index] = spaxel_weight[icube_index] + weight_distance
            spaxel_iflux[icube_index] = spaxel_iflux[icube_index] + 1

#    print('Number of pixels not in ifu cube too low wavelength', ilow)
#    print('Number of pixels not in ifu cube too high wavelength', ihigh)
#    print('Number of pixels not in ifu cube not match', imatch)
# _______________________________________________________________________


def match_det2cube_miripsf(alpha_resol, beta_resol, wave_resol,
                           naxis1, naxis2, naxis3,
                           xcenters, ycenters, zcoord,
                           spaxel_flux,
                           spaxel_weight,
                           spaxel_iflux,
                           spaxel_alpha, spaxel_beta, spaxel_wave,
                           flux,
                           coord1, coord2, wave, alpha_det, beta_det,
                           rois_pixel, roiw_pixel, weight_pixel,
                           softrad_pixel):
    """ Map the detector pixels to the cube spaxels using miri PSF weighting

    Map coordinates coord1,coord2, and wave of the point cloud to which
    spaxels they overlap with in the ifucube. For each spaxel the coord1,coord2
    and wave point cloud members are weighting according to the miri psf and lsf.
    The weighting function is based on the distance the point cloud member
    and spaxel center in the alph-beta coordinate system. The alpha and beta
    value of each point cloud member  and spaxel center is passed to this
    routine.

    Parameters
    ----------
    alpha_resol : numpy.ndarray
      alpha resolution table
    beta_resol : numpy.ndarray
      beta resolution table
    wave_resol : numpy.ndarray
      wavelength resolution table
    naxis1 : int
       size of the ifucube in 1st axis
    naxis2 : int
       size of the ifucube in 2nd axis
    naxis3 : int
       size of the ifucube in 3rd axis
    xcenter : numpy.ndarray
       spaxel center locations 1st dimensions.
    ycenter : numpy.ndarray
       spaxel center locations 2nd dimensions.
    zcoord : numpy.ndarray
        spaxel center locations in 3rd dimensions
    spaxel_flux : numpy.ndarray
       contains the weighted summed detector fluxes that fall withi the roi
    spaxel_weight : numpy.ndarray
       contains the summed weights assocated with the detector fluxes
    spaxel_iflux : numpy.ndarray
       number of detector pixels falling with roi of spaxel center
    flux : numpy.ndarray
       array of detector fluxes associated with each position in
       coorr1, coord2, wave
    spaxel_alpha : numpy.ndarray
       alpha value of spaxel centers
    spaxel_beta : numpy.ndarray
       beta value of spaxel centers
    spaxel_wave : numpy.ndarray
       alpha value of spaxel centers
    coord1 : numpy.ndarray
       contains the spatial coordinate for 1st dimension for the mapped
       detector pixel
    coord2 : numpy.ndarray
       contains the spatial coordinate for 2nd dimension for the mapped
       detector pixel
    wave : numpy.ndarray
       contains the spectral coordinate  for the mapped detector pixel
    alpha_det : numpy.ndarray
       alpha coordinate of mapped detector pixels
    beta_det :  numpy.ndarray
       beta coordinate of mapped detector pixels
    rois_pixel : flaot
       region of influence size in spatial dimension
    roiw_pixel : float
       region of influence size in spectral dimension
    weight_power : float
       msm weighting parameter
    softrad_pxiel :float
       weighting paramter

    Returns
    -------
    spaxel_flux, spaxel_weight, and spaxel_ifux updated with the information
    from the detector pixels that fall within the roi if the spaxel center.

    """

    nplane = naxis1 * naxis2
    # now loop over the pixel values for this region and find the spaxels
    # that fall within the region of interest.
    nn = coord1.size
    #    print('looping over n points mapping to cloud',nn)
    # _______________________________________________________________________
    for ipt in range(0, nn - 1):
        lower_limit = softrad_pixel[ipt]
        # ___________________________________________________________________
        # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.
        # find the spaxels that fall withing ROI of point cloud defined by
        # coord1,coord2,wave

        xdistance = (xcenters - coord1[ipt])
        ydistance = (ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)

        indexr = np.where(radius <= rois_pixel[ipt])
        indexz = np.where(abs(zcoord - wave[ipt]) <= roiw_pixel[ipt])

        # _______________________________________________________________
        # loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for ir, rr in enumerate(indexr[0]):
                # ________________________________________________________
                # if weight is miripsf -distances determined in alpha-beta
                # coordinate system

                weights = FindNormalizationWeights(wave[ipt],
                                                   wave_resol,
                                                   alpha_resol,
                                                   beta_resol)

                cube_index = istart + rr

                alpha_distance = alpha_det[ipt] - spaxel_alpha[cube_index]
                beta_distance = beta_det[ipt] - spaxel_beta[cube_index]
                wave_distance = abs(wave[ipt] - spaxel_wave[cube_index])

                xn = alpha_distance / weights[0]
                yn = beta_distance / weights[1]
                wn = wave_distance / weights[2]

                # only included the spatial dimensions
                wdistance = (xn * xn + yn * yn + wn * wn)
                weight_distance = np.power(np.sqrt(wdistance), weight_pixel[ipt])
# ________________________________________________________________________________
# We have found the weight_distance based on instrument type

                if weight_distance < lower_limit:
                    weight_distance = lower_limit
                weight_distance = 1.0 / weight_distance

                spaxel_flux[cube_index] = spaxel_flux[cube_index] + weight_distance * flux[ipt]
                spaxel_weight[cube_index] = spaxel_weight[cube_index] + weight_distance
                spaxel_iflux[cube_index] = spaxel_iflux[cube_index] + 1
# _______________________________________________________________________


def FindNormalizationWeights(wavelength,
                             wave_resol,
                             alpha_resol,
                             beta_resol):
    """ Routine used in MIRI PSF weighting to normalize data


    Normalize weighting to each point cloud that contributes to the
    ROI of a spaxel. The normalization of weighting is determined from
    width of PSF as wellas wavelength resolution

    Parameters
    ----------
    wave_resol : numpy.ndarray
      wavelength resolution array
    alpha_resol : numpy.ndarray
      alpha psf resolution array
    beta_resol : numpy.ndarray
      beta psf resolution array

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
