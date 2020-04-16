 # Routines used in Spectral Cube Building
import numpy as np
import math
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
#________________________________________________________________________________
def match_det2cube_msm(naxis1, naxis2, naxis3,
                       cdelt1, cdelt2, cdelt3,
                       rois, roiw,
                       msm_weight_power,
                       xcoord, ycoord, zcoord,
                       spaxel_flux,
                       spaxel_weight,
                       spaxel_iflux,
                       spaxel_var,
                       flux,
                       err,
                       coord1, coord2, wave):


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
    xcoord,ycoord: spaxel center locations in 1st and 2nd dimensions. These values 1 D
    zcoord: spaxel center locations in 3rd dimensions
    spaxel_flux: contains the weighted summed detector fluxes that fall withi the roi
    spaxel_weight:  contains the summed weights assocated with the detector fluxes
    spaxel_iflux: number of detector pixels falling with roi of spaxel center
    spaxel_var: contains the weighted summed detector variances that fall withi the roi
    flux: array of detector fluxes associated with each position in  coorr1, coord2, wave
    err: array of detector errors associated with each position in  coorr1, coord2, wave
    coord1, coord2, wave
    Returns
    -------
    spaxel_flux, spaxel_weight, spaxel_ifux, spaxel_var updated with the information from the
    detector pixels that fall within the roi if the spaxel center.
    """

    nplane = naxis1 * naxis2
    lower_limit = 0.01

# now loop over the pixel values for this region and find the spaxels that fall
# withing the region of interest.
    nn = coord1.size

# transform to ifu cube grid system
    c1 = (coord1 - xcoord[0]) / cdelt1
    c2 = (coord2 - ycoord[0]) / cdelt2
    c3 = (wave - zcoord[0]) / cdelt3

    spatial_box = math.ceil(rois / cdelt1)
    spectral_radius = math.ceil(roiw / cdelt3)

    c1_min = (c1 - spatial_box).astype(int)
    c1_max = (c1 + spatial_box).astype(int)
    c1_min[c1_min < 0] = 0
    c1_max[c1_max >= naxis1] = naxis1 - 1

    c2_min = (c2 - spatial_box).astype(int)
    c2_min[c2_min < 0] = 0
    c2_max = (c2 + spatial_box).astype(int)
    c2_max[c2_max >= naxis2] = naxis2 - 1

    c3_min = (c3 - spectral_radius).astype(int)
    c3_min[c3_min < 0] = 0
    c3_max = (c3 + spectral_radius).astype(int)
    c3_max[c3_max >= naxis3] = naxis3 - 1

#    print('looping over n points mapping to cloud',nn)
#________________________________________________________________________________
    for ipt in range(0, nn - 1):
#________________________________________________________________________________
        # xcoord is 1 d array of center of spaxel in x direction
        # ycoord is 1 d array of center of spaxel in y direction
        # find the spaxels that fall withing ROI of point cloud defined  by
        # coord1,coord2,wave
        ixrange = np.arange(c1_min[ipt], c1_max[ipt] + 1)
        iyrange = np.arange(c2_min[ipt], c2_max[ipt] + 1)
        izrange = np.arange(c3_min[ipt], c3_max[ipt] + 1)

        xdistance = (xcoord[ixrange] - coord1[ipt])
        ydistance = (ycoord[iyrange] - coord2[ipt])
        zdistance = (zcoord[izrange] - wave[ipt])

        # form  2D array of xy plane indices
        ixindex = np.tile(ixrange[np.newaxis].T, [1, iyrange.shape[0]])
        iyindex = np.tile(iyrange, [ixindex.shape[0], 1])
        xyindex = (iyindex * naxis1) + ixindex

        # for 2D array of  x distance  y distance values
        xvalues = np.tile(xdistance[np.newaxis].T, [1, ydistance.shape[0]])
        yvalues = np.tile(ydistance, [xvalues.shape[0], 1])
        radius = np.sqrt(xvalues * xvalues + yvalues * yvalues)

        # pull out values that fall in the roi
        indexr = np.where(radius <= rois)
        indexz = np.where(abs(zdistance) <= roiw)

        # grab the index of values that fall in the roi
        xyindex_good = xyindex[indexr]
        zindex_good = izrange[indexz]

        d1 = np.array(xvalues[indexr]) / cdelt1
        d2 = np.array(yvalues[indexr]) / cdelt2
        d3 = np.array(zdistance[indexz]) / cdelt3
        dxy = (d1 * d1) + (d2 * d2)

        # shape of dxy is #indexr or number of overlaps in spatial plane
        # shape of d3 is #indexz or number of overlaps in spectral plane
        # shape of dxy_matrix & d3_matrix  (#indexr, #indexz)
        # rows = number of overlaps in spatial plane
        # cols = number of overlaps in spectral plane
        dxy_matrix = np.tile(dxy[np.newaxis].T, [1, d3.shape[0]])
        d3_matrix = np.tile(d3 * d3, [dxy_matrix.shape[0], 1])

        wdistance = dxy_matrix + d3_matrix
        weight_distance = np.power(np.sqrt(wdistance), msm_weight_power)
        weight_distance[weight_distance < lower_limit] = lower_limit
        weight_distance = 1.0 / weight_distance
        weight_distance = weight_distance.flatten('F')
        weighted_flux = weight_distance * flux[ipt]
        weighted_var = (weight_distance * err[ipt]) * (weight_distance * err[ipt])

        icube_index = [iz * nplane + ixy for iz in zindex_good for ixy in xyindex_good]
        spaxel_flux[icube_index] = spaxel_flux[icube_index] + weighted_flux
        spaxel_weight[icube_index] = spaxel_weight[icube_index] + weight_distance
        spaxel_iflux[icube_index] = spaxel_iflux[icube_index] + 1
        spaxel_var[icube_index] = spaxel_var[icube_index] + weighted_var
#_______________________________________________________________________
