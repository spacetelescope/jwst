 # Routines used in Spectral Cube Building
import numpy as np
import math
import logging
from . import coord

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
#________________________________________________________________________________
def match_det2cube_msm(naxis1, naxis2, naxis3,
                       cdelt1, cdelt2, cdelt3,
                       zcdelt3,
                       roiw,
                       xcenters, ycenters, zcoord,
                       spaxel_flux,
                       spaxel_weight,
                       spaxel_iflux,
                       flux,
                       coord1, coord2, wave,
                       rois_pixel, roiw_pixel, weight_pixel, softrad_pixel):

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
    xcenter,ycenter: spaxel center locations in 1st and 2nd dimensions. These values are 2 X 2 grid
    spaxel centers.
    zcoord: spaxel center locations in 3rd dimensions
    spaxel_flux: contains the weighted summed detector fluxes that fall withi the roi
    spaxel_weight:  contains the summed weights assocated with the detector fluxes
    spaxel_iflux: number of detector pixels falling with roi of spaxel center
    flux: array of detector fluxes associated with each position in  coorr1, coord2, wave
    coord1, coord2, wave
    Returns
    -------
    spaxel_flux, spaxel_weight, and spaxel_ifux updated with the information from the
    detector pixels that fall within the roi if the spaxel center.
    """

    nplane = naxis1 * naxis2
    
# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn = coord1.size

    ilow = 0 
    ihigh = 0 
    imatch  = 0 
#    print('looping over n points mapping to cloud',nn)
#________________________________________________________________________________
    for ipt in range(0, nn - 1):
#________________________________________________________________________________
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
        # the edge cases are skipped and not included in final IFUcube. 
        if len(indexz[0]) == 0:
            if wave[ipt] < zcoord[0]: 
               ilow = ilow + 1

            elif wave[ipt] > zcoord[-1]: 
               ihigh = ihigh + 1
            else:
                imatch = imatch + 1
#                print(' no z match found ',wave[ipt],roiw_pixel[ipt],roiw)
#                print(zcoord[0:10])
#                print(zcoord[naxis3-11:naxis3])
#                print(zcoord[-1])
#                diff = abs(zcoord[naxis3-11:naxis3] - wave[ipt])
#                print(diff)
#                exit()
        else:
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

    print('Number of pixels not in ifu cube too low wavelength', ilow)
    print('Number of pixels not in ifu cube too high wavelength', ihigh)
    print('Number of pixels not in ifu cube not match', imatch)
#_______________________________________________________________________
def match_det2cube_miripsf(alpha_resol, beta_resol, wave_resol,
                           worldtov23,
                           v2ab_transform,
                           naxis1, naxis2, naxis3,
                           cdelt1, cdelt2, cdelt3,
                           crval1, crval2,
                           rois, roiw,
                           weight_power,
                           xcenters, ycenters, zcoord,
                           spaxel_flux,
                           spaxel_weight,
                           spaxel_iflux,
                           flux,
                           coord1, coord2, wave, alpha_det, beta_det):
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
    spaxel_flux: contains the weighted summed detector fluxes that fall withi the roi
    spaxel_weight:  contains the summed weights assocated with the detector fluxes
    spaxel_iflux: number of detector pixels falling with roi of spaxel center
    coord1,coord2,wave pixel coordinates mapped to output frame
    alpha_det,beta_det alpha,beta values of pixel coordinates:

    Returns
    -------
    spaxel_flux, spaxel_weight, and spaxel_ifux updated with the information from the
    detector pixels that fall within the roi if the spaxel center.

    """

    nplane = naxis1 * naxis2
    lower_limit = 0.01


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

        xdistance = (xcenters - coord1[ipt])
        ydistance = (ycenters - coord2[ipt])
        radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)
        indexr = np.where(radius <= rois)
        indexz = np.where(abs(zcoord - wave[ipt]) <= roiw)

        zlam = zcoord[indexz]        # z Cube values falling in wavelength roi
        xi_cube = xcenters[indexr]   # x Cube values within radius
        eta_cube = ycenters[indexr]  # y cube values with the radius
#________________________________________________________________________________
# loop over the points in the ROI
        for iz, zz in enumerate(indexz[0]):
            istart = zz * nplane
            for ir, rr in enumerate(indexr[0]):
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
                spaxel_flux[cube_index] = spaxel_flux[cube_index] + weight_distance * flux[ipt]
                spaxel_weight[cube_index] = spaxel_weight[cube_index] + weight_distance
                spaxel_iflux[cube_index] = spaxel_iflux[cube_index] + 1
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
    print('wavelength,wavecenter',wavelength,wavecenter)
    wave_diff = wavelength - wavecenter
    resolution = a_ave + b_ave * wave_diff + c_ave * wave_diff * wave_diff
    lambda_weight = wavelength / resolution
    weight = [alpha_weight, beta_weight, lambda_weight]
    return weight
