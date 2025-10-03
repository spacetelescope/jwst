import logging
from collections import OrderedDict

import astropy.units as u
import numpy as np
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.table import QTable
from astropy.time import Time, TimeDelta
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture
from photutils.centroids import centroid_sources
from photutils.psf import GaussianPRF, PSFPhotometry
from photutils.utils import CutoutImage
from stdatamodels.jwst.datamodels import CubeModel

__all__ = ["convert_data_units", "tso_aperture_photometry", "tso_source_centroid"]

log = logging.getLogger(__name__)


def convert_data_units(datamodel, gain_2d=None):
    """
    Convert flux units for the input model.

    The datamodel is updated in place.

    Parameters
    ----------
    datamodel : `CubeModel`
        The input `CubeModel` of a TSO imaging observation.
    gain_2d : ndarray or None, optional
        The gain for all pixels.  Required if the input units are "DN/s".
    """
    if datamodel.meta.bunit_data == "MJy/sr":
        # Convert the input data and errors from MJy/sr to Jy
        factor = 1e6 * datamodel.meta.photometry.pixelarea_steradians
        datamodel.data *= factor
        datamodel.err *= factor
        datamodel.meta.bunit_data = "Jy"
        datamodel.meta.bunit_err = "Jy"
    elif datamodel.meta.bunit_data == "DN/s" and gain_2d is not None:
        # Convert the input data and errors from DN/s to electrons
        factor = datamodel.meta.exposure.integration_time * gain_2d
        datamodel.data *= factor
        datamodel.err *= factor
        datamodel.meta.bunit_data = "electron"
        datamodel.meta.bunit_err = "electron"
    else:
        # Unexpected units - leave them as-is
        log.warning(
            f"Unexpected data units: {datamodel.meta.bunit_data}. "
            "Photometry will be produced using the input units."
        )


def tso_aperture_photometry(
    datamodel,
    xcenter,
    ycenter,
    radius,
    radius_inner,
    radius_outer,
    centroid_x=None,
    centroid_y=None,
    psf_width_x=None,
    psf_width_y=None,
    psf_flux=None,
):
    """
    Create a photometric catalog for TSO imaging observations.

    Parameters
    ----------
    datamodel : `CubeModel`
        The input `CubeModel` of a TSO imaging observation.
    xcenter, ycenter : float or ndarray
        The ``x`` and ``y`` center of the aperture.  If a single value
        is provided, it will be used for all integrations.  If an array
        is provided, its size must match the number of integrations in the
        datamodel.
    radius : float
        The radius (in pixels) of the circular aperture.
    radius_inner, radius_outer : float
        The inner and outer radii (in pixels) of the circular-annulus
        aperture, used for local background estimation.
    centroid_x, centroid_y : ndarray or None, optional
        An array of fit centroid values in the x- and y-direction,
        one for each integration. If provided, the arrays will be added
        to the output catalog in the ``centroid_x`` and ``centroid_y``
        columns.
    psf_width_x, psf_width_y : ndarray or None, optional
        An array of fit PSF width values (1-sigma) in the x- and y-direction,
        one for each integration. If provided, the arrays will be added
        to the output catalog in the ``psf_width_x`` and ``psf_width_y``
        columns.
    psf_flux : ndarray or None, optional
        An array of PSF flux values, derived from a Gaussian model,
        one value for each integration. If provided, the array will be added
        to the output catalog in the ``psf_flux`` column.

    Returns
    -------
    catalog : `~astropy.table.QTable`
        An astropy QTable (Quantity Table) containing the source
        photometry.
    """
    if not isinstance(datamodel, CubeModel):
        raise TypeError("The input data model must be a CubeModel.")

    # For the SUB64P subarray with the WLP8 pupil, the circular aperture
    # extends beyond the image and the circular annulus does not have any
    # overlap with the image.  In that case, we simply sum all values
    # in the array and skip the background subtraction.
    sub64p_wlp8 = False
    if datamodel.meta.instrument.pupil == "WLP8" and datamodel.meta.subarray.name == "SUB64P":
        sub64p_wlp8 = True

    # Check for a moving center
    nimg = datamodel.data.shape[0]
    xcenter = np.full(nimg, xcenter) if np.isscalar(xcenter) else xcenter
    ycenter = np.full(nimg, ycenter) if np.isscalar(ycenter) else ycenter

    aperture_sum = []
    aperture_sum_err = []
    annulus_sum = []
    annulus_sum_err = []

    if sub64p_wlp8:
        info = (
            "Photometry measured as the sum of all values in the "
            "subarray.  No background subtraction was performed."
        )

        for i in np.arange(nimg):
            aperture_sum.append(np.nansum(datamodel.data[i, :, :]))
            aperture_sum_err.append(np.sqrt(np.nansum(datamodel.err[i, :, :] ** 2)))
    else:
        info = (
            f"Photometry measured in a circular aperture of r={radius} "
            "pixels.  Background calculated as the mean in a "
            f"circular annulus with r_inner={radius_inner} pixels and "
            f"r_outer={radius_outer} pixels."
        )

        for i in np.arange(nimg):
            phot_aper = CircularAperture((xcenter[i], ycenter[i]), r=radius)
            bkg_aper = CircularAnnulus(
                (xcenter[i], ycenter[i]), r_in=radius_inner, r_out=radius_outer
            )

            aperstats = ApertureStats(
                datamodel.data[i, :, :], phot_aper, error=datamodel.err[i, :, :]
            )

            annstats = ApertureStats(
                datamodel.data[i, :, :], bkg_aper, error=datamodel.err[i, :, :]
            )

            aperture_sum.append(aperstats.sum)
            aperture_sum_err.append(aperstats.sum_err)
            annulus_sum.append(annstats.sum)
            annulus_sum_err.append(annstats.sum_err)

    aperture_sum = np.array(aperture_sum)
    aperture_sum_err = np.array(aperture_sum_err)
    annulus_sum = np.array(annulus_sum)
    annulus_sum_err = np.array(annulus_sum_err)

    # construct metadata for output table
    meta = OrderedDict()
    meta["instrument"] = datamodel.meta.instrument.name
    meta["detector"] = datamodel.meta.instrument.detector
    meta["channel"] = datamodel.meta.instrument.channel
    meta["subarray"] = datamodel.meta.subarray.name
    meta["filter"] = datamodel.meta.instrument.filter
    meta["pupil"] = datamodel.meta.instrument.pupil
    meta["target_name"] = datamodel.meta.target.catalog_name
    meta["apertures"] = info

    # initialize the output table
    tbl = QTable(meta=meta)

    # check for the INT_TIMES table extension
    if datamodel.hasattr("int_times") and datamodel.int_times is not None:
        nrows = len(datamodel.int_times)
    else:
        nrows = 0
        log.warning("The INT_TIMES table in the input file is missing or empty.")

    # load the INT_TIMES table data
    if nrows > 0:
        shape = datamodel.data.shape
        num_integ = 1
        if len(shape) > 2:
            num_integ = shape[0]
        int_start = datamodel.meta.exposure.integration_start
        if int_start is None:
            int_start = 1
            log.warning(f"INTSTART not found; assuming a value of {int_start}")

        # Columns of integration numbers & times of integration from the
        # INT_TIMES table.
        int_num = datamodel.int_times["integration_number"]
        mid_utc = datamodel.int_times["int_mid_MJD_UTC"]
        offset = int_start - int_num[0]  # both are one-indexed
        if offset < 0:
            log.warning(
                "Range of integration numbers in science data extends "
                "outside the range in INT_TIMES table."
            )
            log.warning("Can't use INT_TIMES table.")
            del int_num, mid_utc
            nrows = 0  # flag as bad
        else:
            log.debug("Times are from the INT_TIMES table")
            time_arr = mid_utc[offset : offset + num_integ]
            int_times = Time(time_arr, format="mjd", scale="utc")

    if nrows == 0:
        # No int_times available.
        # Compute integration time stamps on the fly
        log.debug("Times computed from EXPSTART and EFFINTTM")
        dt = datamodel.meta.exposure.integration_time
        n_dt = (
            datamodel.meta.exposure.integration_end - datamodel.meta.exposure.integration_start + 1
        )
        dt_arr = np.arange(1, 1 + n_dt) * dt - (dt / 2.0)
        int_dt = TimeDelta(dt_arr, format="sec")
        int_times = Time(datamodel.meta.exposure.start_time, format="mjd") + int_dt

    # populate table columns
    unit = u.Unit(datamodel.meta.bunit_data)
    tbl["MJD"] = int_times.mjd
    tbl["aperture_sum"] = aperture_sum << unit
    tbl["aperture_sum_err"] = aperture_sum_err << unit

    if not sub64p_wlp8:
        tbl["annulus_sum"] = annulus_sum << unit
        tbl["annulus_sum_err"] = annulus_sum_err << unit

        annulus_mean = annulus_sum / bkg_aper.area
        annulus_mean_err = annulus_sum_err / bkg_aper.area
        aperture_bkg = annulus_mean * phot_aper.area
        aperture_bkg_err = annulus_mean_err * phot_aper.area

        tbl["annulus_mean"] = annulus_mean << unit
        tbl["annulus_mean_err"] = annulus_mean_err << unit

        tbl["aperture_bkg"] = aperture_bkg << unit
        tbl["aperture_bkg_err"] = aperture_bkg_err << unit

        net_aperture_sum = aperture_sum - aperture_bkg
        net_aperture_sum_err = np.sqrt(aperture_sum_err**2 + aperture_bkg_err**2)
        tbl["net_aperture_sum"] = net_aperture_sum << unit
        tbl["net_aperture_sum_err"] = net_aperture_sum_err << unit
    else:
        colnames = [
            "annulus_sum",
            "annulus_sum_err",
            "annulus_mean",
            "annulus_mean_err",
            "aperture_bkg",
            "aperture_bkg_err",
        ]
        for col in colnames:
            tbl[col] = np.full(nimg, np.nan)

        tbl["net_aperture_sum"] = aperture_sum << unit
        tbl["net_aperture_sum_err"] = aperture_sum_err << unit

    # Record aperture center and centroid and PSF fit info if available
    pixel_unit = u.Unit("pix")
    deg_unit = u.Unit("deg")
    tbl["aperture_x"] = xcenter << pixel_unit
    tbl["aperture_y"] = ycenter << pixel_unit

    ra_icrs, dec_icrs = datamodel.meta.wcs(xcenter, ycenter)
    tbl["aperture_ra_icrs"] = ra_icrs << deg_unit
    tbl["aperture_dec_icrs"] = dec_icrs << deg_unit

    if centroid_x is not None and centroid_y is not None:
        tbl["centroid_x"] = centroid_x << pixel_unit
        tbl["centroid_y"] = centroid_y << pixel_unit
    if psf_width_x is not None and psf_width_y is not None:
        tbl["psf_width_x"] = psf_width_x << pixel_unit
        tbl["psf_width_y"] = psf_width_y << pixel_unit
    if psf_flux is not None:
        tbl["psf_flux"] = psf_flux << unit

    return tbl


def _fit_source(data, mask, source_mask, xcenter, ycenter, box_size, fit_psf=False):
    """
    Fit the source in all integrations.

    Parameters
    ----------
    data : ndarray of float
        3D data cube (nimage, ny, nx).
    mask : ndarray of bool
        Mask matching the input data (nimage, ny, nx).  True indicates an invalid pixel.
    source_mask : ndarray of bool
        2D mask for likely source position (ny, nx).
    xcenter, ycenter : float
        Starting guess for the source position.
    box_size : int
        Subimage size to fit.
    fit_psf : bool, optional
        If True and a centroid is successfully fit, the source will be fit
        with a Gaussian model and the PSF width and flux will be returned.

    Returns
    -------
    fit_results : tuple
        If fit_psf is False, a 2-tuple is returned, containing:

        centroid_x : ndarray
            The x center of the source for each integration, zero-indexed.
        centroid_y : ndarray
            The y center of the source for each integration, zero-indexed.

        If fit_psf is True, a 5-tuple is returned, additionally containing:

        psf_width_x : ndarray
            An array of PSF fit widths in the x-direction, one per integration.
        psf_width_y : ndarray
            An array of PSF fit widths in the y-direction, one per integration.
        psf_flux : ndarray
            An array of PSF fit flux, one per integration.
    """
    # Set up output arrays
    nimg = data.shape[0]
    centroid_x = np.full(nimg, np.nan)
    centroid_y = np.full(nimg, np.nan)
    if fit_psf:
        psf_width_x = np.full(nimg, np.nan)
        psf_width_y = np.full(nimg, np.nan)
        psf_flux = np.full(nimg, np.nan)
    else:
        psf_width_x, psf_width_y, psf_flux = None, None, None

    for i in range(nimg):
        image = data[i]

        background = np.nanmedian(image[~source_mask])
        background_sub = image - background

        try:
            centroid = centroid_sources(
                background_sub, xcenter, ycenter, mask=mask[i], box_size=box_size
            )
            centroid_x[i] = centroid[0][0]
            centroid_y[i] = centroid[1][0]
        except ValueError as err:
            # Source could not be centroided. Keep NaN in the output array.
            log.debug(f"Centroid failure in image {i}: {str(err)}")

        # Check for centroid out of range
        ny, nx = image.shape
        if centroid_x[i] < 0 or centroid_x[i] >= nx or centroid_y[i] < 0 or centroid_y[i] >= ny:
            log.debug(
                f"Centroid out of range in image {i}: ({centroid_x[i]},{centroid_y[i]}). "
                f"Setting to NaN."
            )
            centroid_x[i] = np.nan
            centroid_y[i] = np.nan

        if not fit_psf or np.isnan(centroid_x[i]) or np.isnan(centroid_y[i]):
            # Skip PSF calculations
            continue

        # Fit to the PSF at the centroid location
        try:
            psf_width_x[i], psf_width_y[i], psf_flux[i] = _psf_fit_gaussian_prf(
                background_sub, mask[i], box_size, centroid_x[i], centroid_y[i]
            )
        except ValueError as err:
            # Source could not be fit. Keep NaN in the output array.
            log.debug(f"PSF fit failure in image {i}: {str(err)}")

    if fit_psf:
        return centroid_x, centroid_y, psf_width_x, psf_width_y, psf_flux
    else:
        return centroid_x, centroid_y


def _psf_fit_gaussian_prf(data, mask, fit_box_width, xcenter, ycenter):
    """
    Fit a source with a GaussianPRF model, using photutils tools.

    Parameters
    ----------
    data : ndarray of float
        Background subtracted image to fit.
    mask : ndarray of bool
        Mask for bad pixels, matching the data shape. True indicates an invalid pixel.
    fit_box_width : int
        Width of the subimage to use for the fit to the source.
    xcenter, ycenter : float
        Centroid position of the source.

    Returns
    -------
    x_width : float
        The x-width of the Gaussian profile (1-sigma).
    y_width : float
        The y-width of the Gaussian profile (1-sigma).
    flux : float
        The integrated flux contained in the Gaussian profile.
    """
    fit_shape = (fit_box_width, fit_box_width)
    cutout = CutoutImage(data, [ycenter, xcenter], fit_shape)
    cutout_mask = CutoutImage(mask, [ycenter, xcenter], fit_shape)

    # Guess FWHM from box width
    x_fwhm = fit_box_width / 2
    y_fwhm = fit_box_width / 2

    # Guess flux level
    flux = np.sum(cutout.data[~cutout_mask.data])

    # Initial parameters for fix
    init_params = QTable()
    init_params["x"] = [xcenter]
    init_params["y"] = [ycenter]
    init_params["flux"] = [flux]
    init_params["x_fwhm"] = [x_fwhm]
    init_params["y_fwhm"] = [y_fwhm]

    # Integrated 2D Gaussian model
    model = GaussianPRF(flux=flux, x_0=xcenter, y_0=ycenter, x_fwhm=x_fwhm, y_fwhm=y_fwhm)
    model.flux.fixed = False
    model.flux.min = 0.0
    model.flux.max = np.inf
    model.x_0.fixed = True
    model.y_0.fixed = True
    model.x_fwhm.min = 0.0
    model.x_fwhm.max = fit_box_width
    model.x_fwhm.fixed = False
    model.y_fwhm.min = 0.0
    model.y_fwhm.max = fit_box_width
    model.y_fwhm.fixed = False

    phot = PSFPhotometry(model, fit_shape)
    try:
        results = phot(data, mask=mask, init_params=init_params)
    except ValueError as err:
        log.debug("PSF fit failure: %s", str(err))
        return np.nan, np.nan, np.nan

    x_width = gaussian_fwhm_to_sigma * results["x_fwhm_fit"][0]
    y_width = gaussian_fwhm_to_sigma * results["y_fwhm_fit"][0]
    flux = results["flux_fit"][0]
    return x_width, y_width, flux


def tso_source_centroid(
    datamodel, xcenter, ycenter, search_box_width=41, fit_box_width=11, source_radius=4.0
):
    """
    Centroid the source and fit a Gaussian PSF to a subimage.

    For each integration, the source centroid is computed as the
    center-of-mass for the subimage.

    The initial fit to the data uses a wider search box centered on the
    planned position to derive an initial guess for the centroid position.
    A secondary fit uses a smaller subimage to compute the final centroid
    position. The subimage in both cases is background-subtracted prior to
    the fit, from a median value of pixels outside an assumed source radius.

    If the fit is successful, a Gaussian PSF is fit to the subimage at the
    centroid location and the PSF width and flux are reported in the output.

    If the fit fails for all integrations in either the initial or the
    secondary pass, the returned centroid values will be all-NaN arrays and
    the PSF values will be None.  If only some integrations fail, or if
    the centroid succeeds but the Gaussian fits fail, individual values
    within the returned arrays will be set to NaN.

    Parameters
    ----------
    datamodel : `CubeModel`
        The input `CubeModel` of a TSO imaging observation.
    xcenter : float
        Initial guess for the x-center of the source.
    ycenter : float
        Initial guess for the y-center of the source.
    search_box_width : int, optional
        Width of the subimage to use for an initial search for the source.
    fit_box_width : int, optional
        Width of the subimage to use for the final fit to the source.
    source_radius : float, optional
        Expected PSF source radius, used to mask the source for approximate
        background estimation.

    Returns
    -------
    centroid_x : ndarray
        The x center of the source for each integration, zero-indexed.
    centroid_y : ndarray
        The y center of the source for each integration, zero-indexed.
    psf_width_x : ndarray or None
        An array of PSF fit widths in the x-direction, one per integration.
    psf_width_y : ndarray or None
        An array of PSF fit widths in the y-direction, one per integration.
    psf_flux : ndarray or None
        An array of PSF fit flux, one per integration.
    """
    # Get data shape
    nimg, ny, nx = datamodel.data.shape
    yidx, xidx = np.mgrid[:ny, :nx]

    # Mask any pixels marked as DO_NOT_USE or that are NaN
    mask = ((datamodel.dq & 1) != 0) | np.isnan(datamodel.data)

    # We'll get a rough background estimate from each full image,
    # excluding pixels likely to be affected by the source
    source_mask = ((xidx - xcenter) ** 2 + (yidx - ycenter) ** 2) < source_radius**2

    # Get initial centroids from planned position
    centroid_x, centroid_y = _fit_source(
        datamodel.data, mask, source_mask, xcenter, ycenter, search_box_width
    )

    # Check if there were any valid fits
    if not np.any(np.isfinite(centroid_x)) or not np.any(np.isfinite(centroid_y)):
        log.warning("Source could not be centroided.")
        return centroid_x, centroid_y, None, None, None

    # Use the median fitted centroid as the new guess
    xcenter = np.nanmedian(centroid_x)
    ycenter = np.nanmedian(centroid_y)
    log.debug(f"New best guess source position: {xcenter},{ycenter}")

    # Re-centroid with a smaller search box around the new best guess,
    # and fit the PSF at the new centroid location.
    source_mask = ((xidx - xcenter) ** 2 + (yidx - ycenter) ** 2) < source_radius**2
    centroid_x, centroid_y, psf_width_x, psf_width_y, psf_flux = _fit_source(
        datamodel.data, mask, source_mask, xcenter, ycenter, fit_box_width, fit_psf=True
    )

    return centroid_x, centroid_y, psf_width_x, psf_width_y, psf_flux
