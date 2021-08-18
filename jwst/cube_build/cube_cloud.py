""" Map the detector pixels to the cube coordinate system.
This is where the weight functions are used.
"""
import numpy as np
import logging
from shapely.geometry import Polygon
import pdb

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def match_det2cube_driz(naxis1, naxis2, naxis3,
                       cdelt1, cdelt2,
                       zcdelt3,
                       xcenters, ycenters, zcoord,
                       spaxel_flux,
                       spaxel_weight,
                       spaxel_iflux,
                       spaxel_var,
                       flux,
                       err,
                        ccoord, wave, dwave):
    
    """ Map the detector pixels to the cube spaxels using the drizzle parameters
    Match the Point Cloud members to the spaxel centers that fall in the drizzle
    droplet.
    Note that this routine does NOT build the cube by looping over spaxels and
    looking for pixels that contribute to those spaxels.  The runtime is significantly
    better to instead loop over pixels, and look for spaxels that they contribute to.
    This way we can just keep a running sum of the weighted fluxes in a 1-d representation
    of the output cube, which can be normalized by the weights at the end.
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
    spaxel_var: numpy.ndarray
       contains the weighted summed variance within the roi
    flux : numpy.ndarray
       array of detector fluxes associated with each position in
       coorr1, coord2, wave
    err: numpy.ndarray
       array of detector errors associated with each position in
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
    spaxel_flux, spaxel_weight, spaxel_ifux, and spaxel_var updated with the information
    from the detector pixels that fall within the roi of the spaxel center.
    """
    nplane = naxis1 * naxis2

    # Corner coordinates
    xcc1, ycc1 = ccoord[0], ccoord[1]
    xcc2, ycc2 = ccoord[2], ccoord[3]
    xcc3, ycc3 = ccoord[4], ccoord[5]
    xcc4, ycc4 = ccoord[6], ccoord[7]

    xleft=xcenters-cdelt1/2.
    xright=xcenters+cdelt1/2.
    ybot=ycenters-cdelt2/2.
    ytop=ycenters+cdelt2/2.

    # now loop over the pixel values for this region and find the spaxels that fall
    # within the region of interest.

    nn = xcc1.size
    
    print('looping over n pixels',nn, xcc1.size)
    # ________________________________________________________________________________
    for ipt in range(0, nn - 1):
        # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
        # cube coordinates.

        # Z region of interest
        zreg = np.abs(dwave[ipt] + np.max(zcdelt3))
        indexz = np.where(abs(zcoord - wave[ipt]) <= zreg)

        # Output polygons are regular grid squares aligned with coordinate frame, which
        # greatly simplifies search for overlap
        indexr = np.where((xleft < np.max([xcc1[ipt],xcc2[ipt],xcc3[ipt],xcc4[ipt]])) \
                          & (xright > np.min([xcc1[ipt],xcc2[ipt],xcc3[ipt],xcc4[ipt]])) \
                          & (ybot < np.max([ycc1[ipt], ycc2[ipt], ycc3[ipt], ycc4[ipt]])) \
                          & (ytop > np.min([ycc1[ipt], ycc2[ipt], ycc3[ipt], ycc4[ipt]])))

        # Find the cube spectral planes that this input point will contribute to
        if len(indexz[0]) > 0:
            # What is the fractional wavelength overlap?
            ptmin, ptmax = wave[ipt] - dwave[ipt]/2. , wave[ipt] + dwave[ipt]/2.
            spxmin, spxmax = zcoord[indexz] - zcdelt3[indexz]/2. , zcoord[indexz] + zcdelt3[indexz]/2.
            zoverlap = np.zeros(len(indexz[0]))
            for ii in range(0,len(indexz[0])):
                zoverlap[ii] = max(max((spxmax[ii]-ptmin), 0) - max((spxmax[ii]-ptmax), 0) \
                                  - max((spxmin[ii]-ptmin), 0), 0)

            # What is the fractional spatial overlap?
            aoverlap = np.zeros(len(indexr[0]))
            poly1 = Polygon([(xcc1[ipt],ycc1[ipt]), (xcc2[ipt],ycc2[ipt]), \
                             (xcc3[ipt],ycc3[ipt]), (xcc4[ipt],ycc4[ipt])])

            for ii in range(0,len(indexr[0])):
                poly2 = Polygon([(xcenters[indexr][ii] - cdelt1/2., ycenters[indexr][ii] - cdelt2/2.), \
                                 (xcenters[indexr][ii] - cdelt1/2., ycenters[indexr][ii] + cdelt2/2.), \
                                 (xcenters[indexr][ii] + cdelt1/2., ycenters[indexr][ii] + cdelt2/2.), \
                                 (xcenters[indexr][ii] + cdelt1/2., ycenters[indexr][ii] - cdelt2/2.)])
                intersect = poly1.intersection(poly2)
                aoverlap[ii] = intersect.area

            # shape of aoverlap is #indexr or number of overlaps in spatial plane
            # shape of zoverlap is #indexz or number of overlaps in spectral plane
            # shape of a_matrix & z_matrix  (#indexr, #indexz)
            # rows = number of overlaps in spatial plane
            # cols = number of overlaps in spectral plane
            a_matrix = np.tile(aoverlap[np.newaxis].T, [1, zoverlap.shape[0]])
            z_matrix = np.tile(zoverlap, [a_matrix.shape[0], 1])

            weight = a_matrix * z_matrix

            weight = weight.flatten('F')
            weighted_flux = weight * flux[ipt]
            weighted_var = (weight * err[ipt]) * (weight * err[ipt])

            # Identify all of the cube spaxels (ordered in a 1d vector) that this input point contributes to
            icube_index = [iz * nplane + ir for iz in indexz[0] for ir in indexr[0]]

            # Add the weighted flux and variance to running 1d cubes, along with the weights
            # (for later normalization), and point count (for information)
            spaxel_flux[icube_index] = spaxel_flux[icube_index] + weighted_flux
            spaxel_weight[icube_index] = spaxel_weight[icube_index] + weight
            spaxel_var[icube_index] = spaxel_var[icube_index] + weighted_var
            for iw in range(weight.shape[0]):
                if weight[iw] > 0:
                   spaxel_iflux[icube_index[iw]] = spaxel_iflux[icube_index[iw]] + 1

