import logging

import numpy as np
from photutils.centroids import centroid_sources
from astropy.table import QTable
from scipy.optimize import minimize

from stdatamodels.jwst.datamodels import CubeModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def tso_photometry_centroiding(datamodel, gain_model):
    """
    Create a photometric catalog for NIRCam/MIRI TSO imaging observations.

    Parameters
    ----------
    datamodel : `CubeModel`
        The input `CubeModel` of NIRCam/MIRI TSO imaging observation.

    gain_model : `GainModel`
        The gain reference file model.

    Returns
    -------
    centroid_results : numpy.ndarray
        The centroid results from the `TSOCentroidingStep`.
    psfWidth_results : numpy.ndarray
        The PSF width results from the `TSOCentroidingStep`.
    psfFlux_results : numpy.ndarray
        The PSF flux results from the `TSOCentroidingStep`.
    """

    if not isinstance(datamodel, CubeModel):
        raise ValueError('The input data model must be a CubeModel.')

    # Need the FITS WCS X/YREF_SCI values for setting the
    # photometry aperture location
    if datamodel.meta.wcsinfo.siaf_xref_sci is None:
        raise ValueError('XREF_SCI is missing.')
    if datamodel.meta.wcsinfo.siaf_yref_sci is None:
        raise ValueError('YREF_SCI is missing.')

    if datamodel.meta.bunit_data == 'MJy/sr':
        # Convert the input data and errors from MJy/sr to Jy
        factor = 1e6 * datamodel.meta.photometry.pixelarea_steradians
        datamodel.data *= factor
        datamodel.err *= factor
        datamodel.meta.bunit_data = 'Jy'
        datamodel.meta.bunit_err = 'Jy'
    elif datamodel.meta.bunit_data == 'DN/s':
        # Convert the input data and errors from DN/s to electrons
        factor = datamodel.meta.exposure.integration_time*gain_model.data
        datamodel.data *= factor
        datamodel.err *= factor
        datamodel.meta.bunit_data = 'electron'
        datamodel.meta.bunit_err = 'electron'
    elif datamodel.meta.bunit_data == 'DN':
        # Convert the input data and errors from DN to electrons
        factor = gain_model.data
        datamodel.data *= factor
        datamodel.err *= factor
        datamodel.meta.bunit_data = 'electron'
        datamodel.meta.bunit_err = 'electron'
    else:
        # Unexpected units - leave them as-is
        pass

    # Mask any pixels marked as DO_NOT_USE or that are NaN
    mask = datamodel.dq % 2 == 1 | np.isnan(datamodel.data)

    # Get centroid position from FITS header as initial guess
    xcenter_temp = datamodel.meta.wcsinfo.siaf_xref_sci - 1  # 1-based origin
    ycenter_temp = datamodel.meta.wcsinfo.siaf_yref_sci - 1  # 1-based origin

    nimg = datamodel.data.shape[0]

    if datamodel.meta.exposure.type == 'MIR_IMAGE':
        box_radius_broad = 20
        # Start with broad search around an initial guess because the SIAF
        # values in the MIRI FITS header seem to be off for some reason
        centroids = []
        for i in np.arange(nimg):
            # Do an approximate background subtraction to help centroiding
            data_temp = np.copy(datamodel.data[i, :, :])
            mask_temp = np.copy(mask[i, :, :])
            xpixels, ypixels = np.indices(data_temp.shape)
            # Mask pixels likely to be affected by the source to get a less
            # biased background
            maskRadius = 30
            mask_temp[(xpixels-xcenter_temp)**2 + (ypixels-ycenter_temp)**2 < maskRadius] = True
            data_temp = np.ma.masked_where(mask_temp, data_temp)
            data_temp = datamodel.data[i, :, :] - np.ma.median(data_temp)

            centroid = centroid_sources(data_temp,
                                        xcenter_temp, ycenter_temp,
                                        mask=mask[i, :, :], box_size=(box_radius_broad*2+1))
            xcenter_temp, ycenter_temp = np.array(centroid).flatten()
            # Store the outputs
            centroids.append(np.array([xcenter_temp, ycenter_temp]))
        # Convert the list to a numpy array
        centroids = np.array(centroids)

        # Use the median fitted centroid as new guess
        xcenter_temp, ycenter_temp = np.median(centroids, axis=0)

    # Do a centroiding pass with a tight box around the decent centroid guess
    box_radius_precise = 5
    centroids = []
    psfWidths = []
    psfFluxes = []
    for i in np.arange(nimg):
        # Do an approximate background subtraction to help centroiding
        data_temp = np.copy(datamodel.data[i, :, :])
        mask_temp = np.copy(mask[i, :, :])
        xpixels, ypixels = np.indices(data_temp.shape)
        # Mask pixels likely to be affected by the source to get a less
        # biased background
        if datamodel.meta.exposure.type == 'MIR_IMAGE':
            maskRadius = 15
        elif datamodel.meta.exposure.type == 'NRC_TSIMAGE':
            # FINDME: This value is only for defocused photometry. Need to figure out a good value for focused photometry.
            maskRadius = 100
        mask_temp[(xpixels-xcenter_temp)**2 + (ypixels-ycenter_temp)**2 < maskRadius] = True
        data_temp = np.ma.masked_where(mask_temp, data_temp)
        data_centroiding = datamodel.data[i, :, :] - np.ma.median(data_temp)
        data_centroiding = np.ma.masked_where(mask[i, :, :], data_centroiding)

        centroid = centroid_sources(data_centroiding,
                                    xcenter_temp, ycenter_temp,
                                    mask=mask[i, :, :], box_size=(box_radius_precise*2+1))
        x, y = np.array(centroid).flatten()
        # Store the outputs
        centroids.append(np.array([x, y]))

        # Cut out a box around the centroid to measure the PSF width
        minx = -maskRadius+int(np.round(x))
        maxx = maskRadius+int(np.round(x))+1
        miny = -maskRadius+int(np.round(y))
        maxy = maskRadius+int(np.round(y))+1
        data_cutout = np.ma.masked_where(
            mask[i, miny:maxy, minx:maxx],
            datamodel.data[i, miny:maxy, minx:maxx])
        data_cutout -= np.ma.median(data_temp)
        ypixels, xpixels = np.mgrid[miny:maxy, minx:maxx]
        # Estimate the PSF width by calculating sqrt of the second moment
        xWidth = np.ma.sqrt(np.ma.sum((xpixels - x)**2 * data_cutout) /
                            np.ma.sum(data_cutout))/2
        yWidth = np.ma.sqrt(np.ma.sum((ypixels - y)**2 * data_cutout) /
                            np.ma.sum(data_cutout))/2
        y_int, x_int = int(np.round(y)), int(np.round(x))
        amp = np.ma.max(data_centroiding[y_int-1:y_int+2, x_int-1:x_int+2])

        # The initial guess for [Gaussian amplitude, xsigma, ysigma]
        initial_guess = [amp, xWidth, yWidth]
        # The bounds for the fit
        bounds = [(amp/2, amp*2), (xWidth/2, xWidth*2), (yWidth/2, yWidth*2)]

        # Define the function to minimize
        def minfunc(params, data, x_mesh, y_mesh, x, y):
            amp, sx, sy = params
            model = amp*np.ma.exp(-0.5*((x_mesh-x)**2/sx**2 + (y_mesh-y)**2/sy**2))
            return np.ma.mean((data - model)**2)

        # Fit the gaussian width by minimizing minfunc with the Powell method.
        results = minimize(minfunc, initial_guess,
                        args=(data_cutout, xpixels, ypixels, x, y),
                        method='Powell', bounds=bounds)
        if results.success:
            amp, xWidth, yWidth = results.x
        else:
            amp, xWidth, yWidth = np.nan, np.nan, np.nan

        psfWidths.append(np.array([xWidth, yWidth]))
        psfFluxes.append(2*np.pi*amp*xWidth*yWidth)

    # Convert the lists to numpy arrays
    centroids = np.array(centroids)
    psfWidths = np.array(psfWidths)
    psfFluxes = np.array(psfFluxes)

    return centroids, psfWidths, psfFluxes
