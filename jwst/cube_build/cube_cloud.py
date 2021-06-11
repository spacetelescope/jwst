""" Map the detector pixels to the cube coordinate system.
This is where the weight functions are used.
"""
from numba import jit
import numpy as np
import logging
from . import cube_overlap

log = logging.getLogger('numba')
log.setLevel(logging.WARNING)


@jit(nopython=True)
def find_closest_wave(iw,w,
                      wavelength_table,
                      rois_table,
                      roiw_table,
                      softrad_table,
                      weight_power_table,
                      scalerad_table,
                      rois_det,
                      roiw_det,
                      softrad_det,
                      weight_det,
                      scalerad_det):

    """ Given a specific wavelength, find the closest value in the wavelength_table

    """
    ifound = (np.abs(wavelength_table - w)).argmin()
    rois_det[iw] = rois_table[ifound]
    roiw_det[iw] = roiw_table[ifound]
    softrad_det[iw] = softrad_table[ifound]
    weight_det[iw] = weight_power_table[ifound]
    scalerad_det[iw] = scalerad_table[ifound]


@jit(nopython=True)
def map_fov_to_dqplane_miri(start_region, end_region,
                            overlap_partial, overlap_full,
                            naxis1, naxis2,
                            xcenters, ycenters, zcoord,
                            cdelt1, cdelt2,
                            coord1, coord2,
                            wave, roiw_ave, slice_no,
                            spaxel_dq):
    """ Set an initial DQ flag for the IFU cube based on FOV of input data

        Map the FOV of each exposure channel (MIRI) to the DQ plane
        and set an initial DQ flagging. The process is different for MIRI and NIRSpec.
        For MIRI all the slices map roughly to the same FOV across.

        Paramteter
        ---------
        start_region
        end_region
        naxis1
        naxis2
        xcenters
        ycenters
        zcoord
        cdelt1
        cdelt2
        coord1: xi coordinates of input data (~x coordinate in IFU space)
        coord2: eta coordinates of input data (~y coordinate in IFU space)
        wave: wavelength of input data
        roiw_ave: average spectral roi used to determine which wavelength bins
            the input values would be mapped to
        slice_no: integer slice value of input data (used in MIRI case to find
            the points of the edge slices.)

        Returns
        -------
        spaxel_dq
        """

    # MIRI mapping:
    # The FOV is roughly the same for all the wavelength ranges.
    # The offset in the slices makes the calculation of the four corners
    # of the FOV more complicated. So we only use the two slices at
    # the edges of the FOV to define the 4 corners.
    # Note we can not use the wcs.footprint because this footprint only
    # consists of 4 values  ra min, ra max, dec min, dec max and we
    # need 4 corners made of 8 different values.

    # find the wavelength boundaries of the band - use two extreme slices
    wavemin = np.amin(wave)
    wavemax = np.amax(wave)

    # self.zcoord holds the center of the wavelength bin
    imin = (np.abs(zcoord - wavemin)).argmin()
    imax = (np.abs(zcoord - wavemax)).argmin()

    # for each wavelength plane - find the 2 extreme slices to set the FOV
    for w in range(imin, imax):
        wave_distance = np.absolute(zcoord[w] - wave)

        # pull out the two extreme slices in the FOV.
        # use these points to set the FOV: start_region, end_region

        istart = start_region
        coord1_start = None
        index_use = np.where((wave_distance < roiw_ave) & (slice_no == istart))
        if len(index_use[0]) > 0:
            coord2_start = coord2[index_use]
            coord1_start = coord1[index_use]

        iend = end_region
        coord1_end = None
        index_use = np.where((wave_distance < roiw_ave) & (slice_no == iend))
        if len(index_use[0]) > 0:
            coord2_end = coord2[index_use]
            coord1_end = coord1[index_use]

        # if there is valid data on this wavelength plane (not in a gap)
        if coord1_start is not None and coord1_end is not None:

            coord1_total = np.append(np.asarray(coord1_start),np.asarray(coord1_end))
            coord2_total = np.append(np.asarray(coord2_start),np.asarray(coord2_end))

            # from an array of x and y values (contained in coord1_total and coord2_total)
            # determine the footprint
            footprint_all = four_corners(coord1_total, coord2_total)

            isline, footprint = footprint_all
            (xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4) = footprint

            # find the overlap of FOV footprint and with IFU Cube
            xi_corner = np.array([xi1, xi2, xi3, xi4])
            eta_corner = np.array([eta1, eta2, eta3, eta4])

            overlap_fov_with_spaxels(overlap_partial, overlap_full,
                                     cdelt1,cdelt2,
                                     naxis1, naxis2,
                                     xcenters, ycenters,
                                     xi_corner, eta_corner, w, w,
                                     spaxel_dq)


@jit(nopython=True)
def map_fov_to_dqplane_nirspec(overlap_partial,
                               naxis1, naxis2,
                               cdelt1, cdelt2,
                               xcenters, ycenters,
                               xcoord, ycoord, zcoord,
                               coord1, coord2,
                               wave,
                               roiw_ave,
                               slice_no,
                               spaxel_dq):
    """ Set an initial DQ flag for the IFU cube based on FOV of input data

        Map the FOV of each slice (NIRSPEC) to the DQ plane and set an initial DQ
        flagging. For  NIRSpec the 30 different slices map to different FOV on the
        range of wavelengths.

        Paramteter
        ---------
        overlap_partial - intermediate DQ flag
        coord1: xi coordinates of input data (~x coordinate in IFU space)
        coord2: eta coordinates of input data (~y coordinate in IFU space)
        wave: wavelength of input data
        roiw_ave: average spectral roi used to determine which wavelength bins
            the input values would be mapped to
        slice_no: integer slice value of input data (used in MIRI case to find
            the points of the edge slices.)
        """
    # NIRSpec Mapping:
    # The FOV of each NIRSpec slice varies across the wavelength range.
    # Each slice is mapped to each IFU wavelength plane
    # The FOV of the slice is really just a line, so instead of using
    # the routines the finds the overlap between a polygon and regular grid-
    # which is used for MIRI - an algorithm that determines the spaxels that
    # the slice line intersects is used instead.

    # for each of the 30 slices - find the projection of this slice
    # onto each of the IFU wavelength planes.
    for islice in range(30):
        index_slice = np.where(slice_no == islice + 1)

        # find the smaller set of wavelengths to search over for this slice
        wavemin = np.amin(wave[index_slice])
        wavemax = np.amax(wave[index_slice])

        # self.zcoord holds the center of the wavelength bin
        imin = (np.abs(zcoord - wavemin)).argmin()
        imax = (np.abs(zcoord - wavemax)).argmin()

        # loop over valid wavelengths for slice and find the projection of the
        # slice on  the wavelength plane

        for w in range(imin, imax):
            wave_distance = np.absolute(zcoord[w] - wave)
            index_use = np.where((wave_distance < roiw_ave) & (slice_no == islice + 1))
            # using only coord values for wavelength slice
            if len(index_use[0]) > 0:
                coord2_use = coord2[index_use]
                coord1_use = coord1[index_use]
                footprint_all = four_corners(coord1_use, coord2_use)
                isline, footprint = footprint_all
                # isline is true if slice values fall on line (which they will for NIRSpec)
                (xi1, eta1, xi2, eta2, xi3, eta3, xi4, eta4) = footprint
                # find the overlap with IFU Cube
                xi_corner = np.array([xi1, xi2, xi3, xi4])
                eta_corner = np.array([eta1, eta2, eta3, eta4])

                # if isline:
                overlap_slice_with_spaxels(overlap_partial,
                                           cdelt1, cdelt2,
                                           naxis1,
                                           xcoord, ycoord,
                                           xi_corner, eta_corner,
                                           w,
                                           spaxel_dq)
                # else:
                #    overlap_fov_with_spaxels(cdelt1, cdelt2,
                #                             naxis1, naxis2,
                #                             xcenters, ycenters,
                #                             xi_corner, eta_corner, w, w,
                #                             spaxel_dq)
# ********************************************************************************


@jit(nopython=True)
def four_corners(coord1, coord2):
    """ helper function to compute the four corners of the FOV

        From an array of x and y values find the 4 corners enclosing these points
        This routine defines the four corners as
        corner 1: location of min coord2
        corner 2: location of min coord1
        corner 3: location of max coord2
        corner 4: location of max coord1

        Parameter
        ----------
        coord1: array of 4 x corners
        coord2: array of 4 y corners

        Returns
        -------
        Footprint of 4 corners
        Is the data contained in coor1 and coord2 represented by a line if yes:
          isline = True if not isline = False

        """

    isline = False
    index = np.where(coord2 == np.amin(coord2))
    xi_corner1 = coord1[index[0]]
    eta_corner1 = coord2[index[0]]

    index = np.where(coord1 == np.amax(coord1))
    xi_corner2 = coord1[index[0]]
    eta_corner2 = coord2[index[0]]

    index = np.where(coord2 == np.amax(coord2))
    xi_corner3 = coord1[index[0]]
    eta_corner3 = coord2[index[0]]

    index = np.where(coord1 == np.amin(coord1))
    xi_corner4 = coord1[index[0]]
    eta_corner4 = coord2[index[0]]
    footprint = (xi_corner1[0], eta_corner1[0],
                 xi_corner2[0], eta_corner2[0],
                 xi_corner3[0], eta_corner3[0],
                 xi_corner4[0], eta_corner4[0])

    distance_min_points = np.sqrt((xi_corner1 - xi_corner4)**2 +
                                  (eta_corner1 - eta_corner4)**2)

    distance_max_points = np.sqrt((xi_corner2 - xi_corner3)**2 +
                                  (eta_corner2 - eta_corner3)**2)
    dist_tolerance = 0.0001  # tolerance used if points fall on a line
    if ((distance_min_points < dist_tolerance) and (distance_max_points < dist_tolerance)):
        isline = True
    footprint_all = (isline, footprint)

    return footprint_all


@jit(nopython=True)
def overlap_fov_with_spaxels(overlap_partial, overlap_full,
                             cdelt1, cdelt2,
                             naxis1, naxis2,
                             xcenters, ycenters,
                             xi_corner, eta_corner, wmin, wmax,
                             spaxel_dq):

    """find the amount of overlap of FOV with each spaxel

        Given the corners of the FOV  find the spaxels that
        overlap with this FOV.  Set the intermediate spaxel  to
        a value based on the overlap between the FOV for each exposure
        and the spaxel area. The values assigned are:
        a. overlap_partial = overlap partial
        b  overlap_full = overlap_full
        bit_wise combination of these values is allowed to account for
        dithered FOVs.

        Parameter
        ----------
        xi_corner: xi coordinates of the 4 corners of the FOV on the wavelenghth plane
        eta_corner: eta coordinates of the 4 corners of the FOV on the wavelength plane
        wmin: minimum wavelength bin in the IFU cube that this data covers
        wmax: maximum wavelength bin in the IFU cube that this data covers

        Sets
        -------
        spaxel_dq : numpy.ndarray containing intermediate dq flag

        """

    # loop over spaxels in the wavelength plane and set slice_dq
    # lets roughly find the spaxels that might be overlapped

    ximin = np.amin(xi_corner)
    ximax = np.amax(xi_corner)
    etamin = np.amin(eta_corner)
    etamax = np.amax(eta_corner)
    index = np.where(((xcenters - cdelt1 / 2) > ximin) &
                     ((xcenters + cdelt1 / 2) < ximax) &
                     ((ycenters - cdelt2 / 2) > etamin) &
                     ((ycenters + cdelt2 / 2) < etamax))

    wave_slice_dq = np.zeros(naxis2 * naxis1, dtype=np.int32)
    area_box = cdelt1 * cdelt2
    for ixy in range(len(index[0])):
        area_overlap = cube_overlap.sh_find_overlap(xcenters[index][ixy],
                                                    ycenters[index][ixy],
                                                    cdelt1, cdelt2,
                                                    xi_corner, eta_corner)

        overlap_coverage = area_overlap / area_box

        tolerance_dq_overlap = 0.05  # spaxel has to have 5% overlap to flag in FOV

        if overlap_coverage > tolerance_dq_overlap:
            if overlap_coverage > 0.95:
                wave_slice_dq[ixy] = overlap_full
            else:
                wave_slice_dq[ixy] = overlap_partial

    # set for a range of wavelengths
    if wmin != wmax:
        spaxel_dq[wmin:wmax, :] = np.bitwise_or(spaxel_dq[wmin:wmax, :],
                                                wave_slice_dq)

        # set for a single wavelength
    else:
        spaxel_dq[wmin, :] = np.bitwise_or(spaxel_dq[wmin, :],
                                           wave_slice_dq)


@jit(nopython=True)
def overlap_slice_with_spaxels(overlap_partial,
                               cdelt1, cdelt2,
                               naxis1,
                               xcoord, ycoord,
                               xi_corner, eta_corner, w,
                               spaxel_dq):
    """ Set the initial dq plane of indicating if the input data falls on a spaxel

        This algorithm assumes the input data falls on a line in the IFU cube, which is
        the case for NIRSpec slices. The NIRSpec slice's endpoints are used to determine
        which IFU spaxels the slice falls on to set an initial dq flag.

        Parameters
        ---------
        overlap_partial: intermediate dq flag
        xi_corner: holds the x starting and ending points of the slice
        eta_corner: holds the y starting and ending points of the slice
        wavelength: the wavelength bin of the IFU cube working with

        Sets
        ----
        self.spaxel_dq : numpy.ndarray containing intermediate dq flag

        """

    points = findpoints_on_slice(cdelt1, cdelt2,
                                 xcoord, ycoord,
                                 xi_corner, eta_corner)
    num = len(points)

    for i in range(num):
        xpt, ypt = points[i]
        index = (ypt * naxis1) + xpt
        spaxel_dq[w, index] = overlap_partial


@jit(nopython=True)
def findpoints_on_slice(cdelt1, cdelt2,
                        xcoord, ycoord,
                        xi_corner, eta_corner):
    """ Bresenham's Line Algorithm to find points a line intersects with grid.

        Given the endpoints of a line find the spaxels this line intersects.

        Parameters
        -----------
        xi_corner: holds the started in ending x values
        eta_corner: holds the started in ending y values

        Returns
        -------
        Points: a tuple of x,y spaxel values that this line intersects

        """

    # set up line - convert to integer values
    x1 = int((xi_corner[0] - xcoord[0]) / cdelt1)
    y1 = int((eta_corner[0] - ycoord[0]) / cdelt2)
    x2 = int((xi_corner[1] - xcoord[0]) / cdelt1)
    y2 = int((eta_corner[1] - ycoord[0]) / cdelt2)

    dx = x2 - x1
    dy = y2 - y1

    # how steep is it
    is_steep = abs(dy) > abs(dx)

    # Rotate line
    if is_steep:
        x1, y1 = y1, x1
        x2, y2 = y2, x2

    # Swap start and end points if necessary and store swap state
    swapped = False
    if x1 > x2:
        x1, x2 = x2, x1
        y1, y2 = y2, y1
        swapped = True

    # Recalculate differences
    dx = x2 - x1
    dy = y2 - y1

    # calculate error
    error = int(dx / 2.0)
    ystep = -1
    if y1 < y2:
        ystep = 1

    # iterate over grid to generate points between the start and end of line
    y = y1
    points = []
    for x in range(x1, x2 + 1):
        coord = (y, x) if is_steep else (x, y)
        points.append(coord)
        error -= abs(dy)
        if error < 0:
            y += ystep
            error += dx

    # If coords were swapped then reverse
    if swapped:
        points.reverse()

    # iterate over grid to generate points between the start and end of line
    # y = y1
    # points2 = np.array([np.array(list)])
    # for x in range(x1, x2 + 1):
    #     coord = (y, x) if is_steep else (x, y)
    #     points2 = np.append(points2, coord)
    #     error -= abs(dy)
    #     if error < 0:
    #         y += ystep
    #         error += dx

    # points2 = np.reshape(points2, (-1,2))
    # If coords were swapped then reverse
    # if swapped:
    #     points2 = points2[:,:,-1]

    points = np.array(points)
    return points


@jit(nopython=True)
def match_det2cube_emsm(ipt,
                        nplane,
                        cdelt1, cdelt2,
                        zcdelt3,
                        xcenters, ycenters, zcoord,
                        spaxel_flux,
                        spaxel_weight,
                        spaxel_iflux,
                        spaxel_var,
                        flux,
                        err,
                        coord1, coord2, wave,
                        rois_pixel, roiw_pixel,
                        scalerad_pixel):

    """ Map the detector pixels to the cube spaxels using the MSM parameters


    Match the Point Cloud members to the spaxel centers that fall in the ROI.
    For each spaxel the coord1,coord1 and wave point cloud members are weighed
    according to modified shepard method of inverse weighting based on the
    distance between the point cloud member and the spaxel center.

    Note that this routine does NOT build the cube by looping over spaxels and
    looking for pixels that contribute to those spaxels.  The runtime is significantly
    better to instead loop over pixels, and look for spaxels that they contribute to.
    This way we can just keep a running sum of the weighted fluxes in a 1-d representation
    of the output cube, which can be normalized by the weights at the end.

    Parameters
    ----------
    ipt : int
       index of point cloud member
    nplane : int
       naxis1 * naxis 2
    cdelt1 : float
       ifucube spaxel size in axis 1 dimension
    cdelt2 : float
       ifucube spaxel size in axis 2 dimension
    zcdelt3 :float
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
    rois_pixel : float
       region of influence size in spatial dimension
    roiw_pixel : float
       region of influence size in spectral dimension

    Returns
    -------
    spaxel_flux, spaxel_weight, spaxel_ifux, and spaxel_var updated with the information
    from the detector pixels that fall within the roi of the spaxel center.
    """

    # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
    # cube coordinates.
    # find the spaxels that fall withing ROI of point cloud defined  by
    # coord1,coord2,wave

    xdistance = (xcenters - coord1[ipt])
    ydistance = (ycenters - coord2[ipt])
    radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)

    indexr = np.where(radius <= rois_pixel[ipt])
    indexz = np.where(np.absolute(zcoord - wave[ipt]) <= roiw_pixel[ipt])

    # Find the cube spectral planes that this input point will contribute to
    if len(indexz[0]) > 0:

        d1 = (coord1[ipt] - xcenters[indexr]) / cdelt1
        d2 = (coord2[ipt] - ycenters[indexr]) / cdelt2
        d3 = (wave[ipt] - zcoord[indexz]) / zcdelt3[indexz]

        dxy = (d1 * d1) + (d2 * d2)

        # shape of dxy is #indexr or number of overlaps in spatial plane
        # shape of d3 is #indexz or number of overlaps in spectral plane
        # shape of dxy_matrix & d3_matrix  (#indexr, #indexz)
        # rows = number of overlaps in spatial plane
        # cols = number of overlaps in spectral plane

        d32 = d3 * d3
        d3_matrix2 = np.zeros((dxy.shape[0], d3.shape[0]))
        dxy_matrix2 = np.zeros((dxy.shape[0], d3.shape[0]))

        dxy_matrix2 = np.repeat(dxy, d3.shape[0])
        dxy_matrix2 = np.reshape(dxy_matrix2, (dxy.shape[0], d3.shape[0]))

        # we cannot form the d3_matrix using np.repeat because we would need
        # to do the following steps and np.repeat (with axis option is not
        # allowed in numba
        # d33 = np.reshape(d32, (-1, d3.shape[0]))
        # d3_matrix2 =np.repeat(d33,dxy.shape[0],axis=0)

        d3_matrix2[:,0:d3.shape[0]] = d32
        wdistance = dxy_matrix2 + d3_matrix2

        weight_distance = np.exp(-wdistance / (scalerad_pixel[ipt] / cdelt1))
        weight_distance2 = np.transpose(weight_distance)
        weight_distance2 = weight_distance2.flatten()

        weighted_flux = weight_distance2 * flux[ipt]
        weighted_var = (weight_distance2 * err[ipt]) * (weight_distance2 * err[ipt])

        # Identify all of the cube spaxels (ordered in a 1d vector) that this input point contributes to
        icube_index = [iz * nplane + ir for iz in indexz[0] for ir in indexr[0]]

        icube_index = np.array(icube_index)

        # Add the weighted flux and variance to running 1d cubes, along with the weights
        # (for later normalization), and point count (for information)
        spaxel_flux[icube_index] = spaxel_flux[icube_index] + weighted_flux
        spaxel_weight[icube_index] = spaxel_weight[icube_index] + weight_distance2
        spaxel_iflux[icube_index] = spaxel_iflux[icube_index] + 1
        spaxel_var[icube_index] = spaxel_var[icube_index] + weighted_var


@jit(nopython=True)
def match_det2cube_msm(ipt,
                       nplane,
                       cdelt1, cdelt2,
                       zcdelt3,
                       xcenters, ycenters, zcoord,
                       spaxel_flux,
                       spaxel_weight,
                       spaxel_iflux,
                       spaxel_var,
                       flux,
                       err,
                       coord1, coord2, wave,
                       weight_pixel,
                       rois_pixel, roiw_pixel,
                       softrad_pixel):

    """ Map the detector pixels to the cube spaxels using the MSM parameters


    Match the Point Cloud members to the spaxel centers that fall in the ROI.
    For each spaxel the coord1,coord1 and wave point cloud members are weighed
    according to modified shepard method of inverse weighting based on the
    distance between the point cloud member and the spaxel center.

    Note that this routine does NOT build the cube by looping over spaxels and
    looking for pixels that contribute to those spaxels.  The runtime is significantly
    better to instead loop over pixels, and look for spaxels that they contribute to.
    This way we can just keep a running sum of the weighted fluxes in a 1-d representation
    of the output cube, which can be normalized by the weights at the end.

    Parameters
    ----------
    ipt : int
       index of point cloud element
    nplane : int
       naxis1 * naxis 2
    cdelt1 : float
       ifucube spaxel size in axis 1 dimension
    cdelt2 : float
       ifucube spaxel size in axis 2 dimension
    zcdelt3 :float
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
    rois_pixel : float
       region of influence size in spatial dimension
    roiw_pixel : float
       region of influence size in spectral dimension

    Returns
    -------
    spaxel_flux, spaxel_weight, spaxel_ifux, and spaxel_var updated with the information
    from the detector pixels that fall within the roi of the spaxel center.
    """

    # xcenters, ycenters is a flattened 1-D array of the 2 X 2 xy plane
    # cube coordinates.
    # find the spaxels that fall withing ROI of point cloud defined  by
    # coord1,coord2,wave
    lower_limit = softrad_pixel[ipt]
    xdistance = (xcenters - coord1[ipt])
    ydistance = (ycenters - coord2[ipt])
    radius = np.sqrt(xdistance * xdistance + ydistance * ydistance)

    indexr = np.where(radius <= rois_pixel[ipt])
    indexz = np.where(np.absolute(zcoord - wave[ipt]) <= roiw_pixel[ipt])

    # Find the cube spectral planes that this input point will contribute to
    if len(indexz[0]) > 0:

        d1 = (coord1[ipt] - xcenters[indexr]) / cdelt1
        d2 = (coord2[ipt] - ycenters[indexr]) / cdelt2
        d3 = (wave[ipt] - zcoord[indexz]) / zcdelt3[indexz]

        dxy = (d1 * d1) + (d2 * d2)

        # shape of dxy is #indexr or number of overlaps in spatial plane
        # shape of d3 is #indexz or number of overlaps in spectral plane
        # shape of dxy_matrix & d3_matrix  (#indexr, #indexz)
        # rows = number of overlaps in spatial plane
        # cols = number of overlaps in spectral plane

        d32 = d3 * d3
        d3_matrix2 = np.zeros((dxy.shape[0], d3.shape[0]))
        dxy_matrix2 = np.zeros((dxy.shape[0], d3.shape[0]))

        dxy_matrix2 = np.repeat(dxy, d3.shape[0])
        dxy_matrix2 = np.reshape(dxy_matrix2, (dxy.shape[0], d3.shape[0]))

        # we cannot form the d3_matrix using np.repeat because we would need
        # to do the following steps and np.repeat (with axis option is not
        # allowed in numba
        # d33 = np.reshape(d32, (-1, d3.shape[0]))
        # d3_matrix2 =np.repeat(d33,dxy.shape[0],axis=0)

        d3_matrix2[:,0:d3.shape[0]] = d32
        wdistance = dxy_matrix2 + d3_matrix2
        weight_distance22 = np.power(np.sqrt(wdistance), weight_pixel[ipt])
        weight_distance22 = np.transpose(weight_distance22)
        weight_distance22 = weight_distance22.flatten()
        weight_distance22[weight_distance22 < lower_limit] = lower_limit
        weight_distance2 = 1.0 / weight_distance22

        weighted_flux = weight_distance2 * flux[ipt]
        weighted_var = (weight_distance2 * err[ipt]) * (weight_distance2 * err[ipt])

        # Identify all of the cube spaxels (ordered in a 1d vector) that this input point contributes to
        icube_index = [iz * nplane + ir for iz in indexz[0] for ir in indexr[0]]

        icube_index = np.array(icube_index)

        # Add the weighted flux and variance to running 1d cubes, along with the weights
        # (for later normalization), and point count (for information)
        spaxel_flux[icube_index] = spaxel_flux[icube_index] + weighted_flux
        spaxel_weight[icube_index] = spaxel_weight[icube_index] + weight_distance2
        spaxel_iflux[icube_index] = spaxel_iflux[icube_index] + 1
        spaxel_var[icube_index] = spaxel_var[icube_index] + weighted_var
