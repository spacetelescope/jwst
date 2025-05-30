import logging
import warnings

import numpy as np
from astropy import stats
from astropy.stats import sigma_clipped_stats as sigclip
from photutils.aperture import (
    CircularAperture,
    CircularAnnulus,
    RectangularAperture,
    aperture_photometry,
)
from photutils.detection import DAOStarFinder
from scipy.interpolate import interp1d

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags

from jwst.assign_wcs.util import compute_scale
from jwst.extract_1d import spec_wcs
from jwst.extract_1d.apply_apcorr import select_apcorr
from jwst.extract_1d.extract import read_apcorr_ref
from jwst.residual_fringe import utils as rfutils

__all__ = ["ifu_extract1d"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# This is intended to be larger than any possible distance (in pixels)
# between the target and any point in the image; used by locn_from_wcs().
HUGE_DIST = 1.0e10


def ifu_extract1d(
    input_model,
    ref_file,
    source_type,
    subtract_background,
    bkg_sigma_clip,
    apcorr_ref_file=None,
    center_xy=None,
    ifu_autocen=False,
    ifu_rfcorr=False,
    ifu_rscale=None,
    ifu_covar_scale=1.0,
):
    """
    Extract a 1D spectrum from an IFU cube.

    Parameters
    ----------
    input_model : JWST data model for an IFU cube (IFUCubeModel)
        The input model.
    ref_file : str
        File name for the extract1d reference file, in ASDF format.
    source_type : str
        "POINT" or "EXTENDED"
    subtract_background : bool or None
        User supplied flag indicating whether the background should be subtracted.
        If None, the value in the extract_1d reference file will be used.
        If not None, this parameter overrides the value in the
        extract_1d reference file.
    bkg_sigma_clip : float
        Background sigma clipping value to use to remove noise/outliers in background
    apcorr_ref_file : str or None, optional
        File name for aperture correction reference file.
    center_xy : float or None, optional
        A list of 2 pixel coordinate values at which to place the center
        of the extraction aperture for IFU data, overriding any centering
        done by the step.  Two values, in x,y order, are used for extraction
        from IFU cubes. Default is None.
    ifu_autocen : bool, optional
        Switch to turn on auto-centering for point source spectral extraction
        in IFU mode.  Default is False.
    ifu_rfcorr : bool, optional
        Switch to select whether or not to apply a 1d residual fringe correction
        for MIRI MRS IFU spectra.  Default is False.
    ifu_rscale : float, optional
        For MRS IFU data, a value for changing the extraction radius. The value
        provided is the number of PSF FWHMs to use for the extraction radius.
        Values accepted are between 0.5 to 3.0. The default extraction size is
        set to 2 * FWHM. Values below 2 will result in a smaller radius, a value
        of 2 results in no change to radius and a value above 2 results in a larger
        extraction radius.
    ifu_covar_scale : float, optional
        Scaling factor by which to multiply the ERR values in extracted spectra to
        account for covariance between adjacent spaxels in the IFU data cube.

    Returns
    -------
    output_model : MultiSpecModel
        This will contain the extracted spectrum.
    """
    if not isinstance(input_model, datamodels.IFUCubeModel):
        log.error("Expected an IFU cube.")
        raise TypeError("Expected an IFU cube.")

    instrument = input_model.meta.instrument.name
    if instrument is not None:
        instrument = instrument.upper()
    if source_type is not None:
        source_type = source_type.upper()
    if source_type != "POINT" and source_type != "EXTENDED":
        default_source_type = "EXTENDED"
        log.warning(f"Source type was '{source_type}'; setting to '{default_source_type}'.")
        source_type = default_source_type
    else:
        log.info(f"Source type = {source_type}")

    if input_model.meta.instrument.name == "MIRI":
        output_model = datamodels.MRSMultiSpecModel()
        spec_dtype = datamodels.MRSSpecModel().spec_table.dtype
    else:
        output_model = datamodels.MultiSpecModel()
        spec_dtype = datamodels.SpecModel().spec_table.dtype

    output_model.update(input_model, only="PRIMARY")

    slitname = input_model.meta.exposure.type
    if slitname is None:
        slitname = "ANY"

    extract_params = get_extract_parameters(ref_file, bkg_sigma_clip)

    # Add info about IFU auto-centroiding, residual fringe correction and extraction radius scale
    # to extract_params for use later
    extract_params["ifu_autocen"] = ifu_autocen
    extract_params["ifu_rfcorr"] = ifu_rfcorr
    extract_params["ifu_rscale"] = ifu_rscale
    extract_params["ifu_covar_scale"] = ifu_covar_scale

    # If the user supplied extraction center coords,
    # load them into extract_params for use later.
    extract_params["x_center"] = None
    extract_params["y_center"] = None

    if center_xy is not None:
        if len(center_xy) == 2:
            extract_params["x_center"] = float(center_xy[0])
            extract_params["y_center"] = float(center_xy[1])
            log.info(f"Using user-supplied x_center={center_xy[0]}, y_center={center_xy[1]}")
        else:
            log.warning("Incorrect number of values in center_xy; should be two.")

    if subtract_background is not None:
        if subtract_background and source_type == "EXTENDED":
            subtract_background = False
            log.info("Turning off background subtraction because the source is extended.")
        extract_params["subtract_background"] = subtract_background

    (
        ra,
        dec,
        wavelength,
        temp_flux,
        f_var_poisson,
        f_var_rnoise,
        f_var_flat,
        background,
        b_var_poisson,
        b_var_rnoise,
        b_var_flat,
        npixels,
        dq,
        npixels_bkg,
        radius_match,
        x_center,
        y_center,
    ) = extract_ifu(input_model, source_type, extract_params)

    npixels_temp = np.where(npixels > 0.0, npixels, 1.0)
    npixels_bkg_temp = np.where(npixels_bkg > 0.0, npixels_bkg, 1.0)

    # Convert the sum to an average, for surface brightness.
    # Remember variance requires dividing by N^2
    surf_bright = temp_flux / npixels_temp
    sb_var_poisson = f_var_poisson / npixels_temp / npixels_temp
    sb_var_rnoise = f_var_rnoise / npixels_temp / npixels_temp
    sb_var_flat = f_var_flat / npixels_temp / npixels_temp
    background = background / npixels_bkg_temp
    b_var_poisson = b_var_poisson / npixels_bkg_temp / npixels_bkg_temp
    b_var_rnoise = b_var_rnoise / npixels_bkg_temp / npixels_bkg_temp
    b_var_flat = b_var_flat / npixels_bkg_temp / npixels_bkg_temp

    del npixels_temp
    del npixels_bkg_temp

    temp_flux_rf = None
    surf_bright_rf = None
    background_rf = None

    # If selected, apply 1d residual fringe correction to the extracted spectrum

    if (input_model.meta.instrument.name == "MIRI") and (extract_params["ifu_rfcorr"] is True):
        log.info("Applying 1d residual fringe correction.")
        # Determine which MRS channel the spectrum is from
        thischannel = input_model.meta.instrument.channel
        # Valid single-channel values
        validch = ["1", "2", "3", "4"]

        # If a valid single channel, specify it in call to residual fringe code
        if thischannel in validch:
            channel = int(thischannel)
        else:
            channel = ""

        # Embed all calls to residual fringe code in try/except blocks as the default behavior
        # if problems are encountered should be to not apply this optional step

        # Apply residual fringe to the flux array
        try:
            temp_flux_rf = rfutils.fit_residual_fringes_1d(
                temp_flux, wavelength, channel=channel, dichroic_only=False, max_amp=None
            )
        except Exception:
            log.info("Flux residual fringe correction failed- skipping.")

        # Apply residual fringe to the surf_bright array
        try:
            surf_bright_rf = rfutils.fit_residual_fringes_1d(
                surf_bright, wavelength, channel=channel, dichroic_only=False, max_amp=None
            )
        except Exception:
            log.info("Surf bright residual fringe correction failed- skipping.")

        # Apply residual fringe to the background array
        try:
            background_rf = rfutils.fit_residual_fringes_1d(
                background, wavelength, channel=channel, dichroic_only=False, max_amp=None
            )
        except Exception:
            log.info("Background residual fringe correction failed- skipping.")

    pixel_solid_angle = input_model.meta.photometry.pixelarea_steradians
    if pixel_solid_angle is None:
        log.warning("Pixel area (solid angle) is not populated; the flux will not be correct.")
        pixel_solid_angle = 1.0

    # Convert flux from MJy / steradian to Jy.
    flux = temp_flux * pixel_solid_angle * 1.0e6
    f_var_poisson *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2
    f_var_rnoise *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2
    f_var_flat *= pixel_solid_angle**2 * 1.0e12  # (MJy / sr)**2 --> Jy**2

    # surf_bright and background were computed above
    del temp_flux

    error = np.sqrt(f_var_poisson + f_var_rnoise + f_var_flat)
    sb_error = np.sqrt(sb_var_poisson + sb_var_rnoise + sb_var_flat)
    berror = np.sqrt(b_var_poisson + b_var_rnoise + b_var_flat)

    # If we only used the Poisson variance array as a vehicle to pass through
    # non-differentiated errors, clear it again here so that only the total
    # errors pass out into the 1d spectra files.
    if not hasattr(input_model, "var_poisson"):
        f_var_poisson[:] = 0
        sb_var_poisson[:] = 0
        b_var_poisson[:] = 0

    # Deal with covariance in the IFU cube by multiplying 1d spectra errors by a scaling factor
    if extract_params["ifu_covar_scale"] != 1.0:
        log.info("Applying ERR covariance prefactor of %g", extract_params["ifu_covar_scale"])
        error *= extract_params["ifu_covar_scale"]
        sb_error *= extract_params["ifu_covar_scale"]
        berror *= extract_params["ifu_covar_scale"]
        f_var_poisson *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        f_var_rnoise *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        f_var_flat *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        sb_var_poisson *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        sb_var_rnoise *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        sb_var_flat *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        b_var_poisson *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        b_var_rnoise *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]
        b_var_flat *= extract_params["ifu_covar_scale"] * extract_params["ifu_covar_scale"]

    if input_model.meta.instrument.name == "MIRI":
        if temp_flux_rf is None:
            flux_rf = np.full_like(flux, np.nan)
        else:
            flux_rf = temp_flux_rf * pixel_solid_angle * 1.0e6
        del temp_flux_rf

        if background_rf is None:
            background_rf = np.full_like(background, np.nan)
        if surf_bright_rf is None:
            surf_bright_rf = np.full_like(surf_bright, np.nan)

        otab = np.array(
            list(
                zip(
                    wavelength,
                    flux,
                    error,
                    f_var_poisson,
                    f_var_rnoise,
                    f_var_flat,
                    surf_bright,
                    sb_error,
                    sb_var_poisson,
                    sb_var_rnoise,
                    sb_var_flat,
                    dq,
                    background,
                    berror,
                    b_var_poisson,
                    b_var_rnoise,
                    b_var_flat,
                    npixels,
                    flux_rf,
                    surf_bright_rf,
                    background_rf,
                    strict=False,
                )
            ),
            dtype=spec_dtype,
        )
        spec = datamodels.MRSSpecModel(spec_table=otab)

    else:  # NIRSPEC
        otab = np.array(
            list(
                zip(
                    wavelength,
                    flux,
                    error,
                    f_var_poisson,
                    f_var_rnoise,
                    f_var_flat,
                    surf_bright,
                    sb_error,
                    sb_var_poisson,
                    sb_var_rnoise,
                    sb_var_flat,
                    dq,
                    background,
                    berror,
                    b_var_poisson,
                    b_var_rnoise,
                    b_var_flat,
                    npixels,
                    strict=False,
                )
            ),
            dtype=spec_dtype,
        )
        spec = datamodels.SpecModel(spec_table=otab)

    spec.meta.wcs = spec_wcs.create_spectral_wcs(ra, dec, wavelength)

    spec.spec_table.columns["wavelength"].unit = "um"
    spec.spec_table.columns["flux"].unit = "Jy"
    spec.spec_table.columns["flux_error"].unit = "Jy"
    spec.spec_table.columns["flux_var_poisson"].unit = "Jy^2"
    spec.spec_table.columns["flux_var_rnoise"].unit = "Jy^2"
    spec.spec_table.columns["flux_var_flat"].unit = "Jy^2"
    spec.spec_table.columns["surf_bright"].unit = "MJy/sr"
    spec.spec_table.columns["sb_error"].unit = "MJy/sr"
    spec.spec_table.columns["sb_var_poisson"].unit = "(MJy/sr)^2"
    spec.spec_table.columns["sb_var_rnoise"].unit = "(MJy/sr)^2"
    spec.spec_table.columns["sb_var_flat"].unit = "(MJy/sr)^2"
    spec.spec_table.columns["background"].unit = "MJy/sr"
    spec.spec_table.columns["bkgd_error"].unit = "MJy/sr"
    spec.spec_table.columns["bkgd_var_poisson"].unit = "(MJy/sr)^2"
    spec.spec_table.columns["bkgd_var_rnoise"].unit = "(MJy/sr)^2"
    spec.spec_table.columns["bkgd_var_flat"].unit = "(MJy/sr)^2"
    if input_model.meta.instrument.name == "MIRI":
        spec.spec_table.columns["rf_flux"].unit = "Jy"
        spec.spec_table.columns["rf_surf_bright"].unit = "MJy/sr"
        spec.spec_table.columns["rf_background"].unit = "MJy/sr"

    spec.slit_ra = ra
    spec.slit_dec = dec

    if slitname is not None and slitname != "ANY":
        spec.name = slitname
    spec.detector = input_model.meta.instrument.detector

    spec.source_type = source_type
    spec.extraction_x = x_center
    spec.extraction_y = y_center

    if source_type == "POINT" and apcorr_ref_file is not None and apcorr_ref_file != "N/A":
        apcorr_ref_model = read_apcorr_ref(apcorr_ref_file, input_model.meta.exposure.type)

        log.info("Applying Aperture correction.")
        if instrument == "NIRSPEC":
            wl = np.median(wavelength)
        else:
            wl = wavelength.min()

        apcorr = select_apcorr(input_model)(input_model, apcorr_ref_model, location=(ra, dec, wl))

        # determine apcor function and apcor radius to use at each wavelength
        apcorr.match_wavelengths(wavelength)

        # at each IFU wavelength we have the extraction radius
        # defined by radius_match (radius size in pixels)
        for i in range(wavelength.size):
            radius = radius_match[i]
            apcorr.find_apcorr_func(i, radius)

        apcorr.apply(spec.spec_table)

    output_model.spec.append(spec)

    # See output_model.spec[0].meta.wcs instead.
    output_model.meta.wcs = None

    return output_model


def get_extract_parameters(ref_file, bkg_sigma_clip):
    """
    Read extraction parameters for an IFU.

    Parameters
    ----------
    ref_file : dict
        File name for the extract1d reference file, in ASDF format
    bkg_sigma_clip : float
        Background sigma clipping value to use to remove noise/outliers in background.

    Returns
    -------
    dict
        The extraction parameters.
    """
    extract_params = {}

    # for consistency put the bkg_sigma_clip in dictionary: extract_params
    extract_params["bkg_sigma_clip"] = bkg_sigma_clip

    refmodel = datamodels.Extract1dIFUModel(ref_file)
    subtract_background = refmodel.meta.subtract_background
    method = refmodel.meta.method
    subpixels = refmodel.meta.subpixels

    data = refmodel.data
    wavelength = data.wavelength
    radius = data.radius
    inner_bkg = data.inner_bkg
    outer_bkg = data.outer_bkg

    extract_params["subtract_background"] = bool(subtract_background)
    extract_params["method"] = method
    extract_params["subpixels"] = subpixels
    extract_params["wavelength"] = wavelength
    extract_params["radius"] = radius
    extract_params["inner_bkg"] = inner_bkg
    extract_params["outer_bkg"] = outer_bkg

    return extract_params


def extract_ifu(input_model, source_type, extract_params):
    """
    Perform 1D extraction for IFU data.

    Parameters
    ----------
    input_model : IFUCubeModel
        The input model.
    source_type : str
        "POINT" or "EXTENDED"
    extract_params : dict
        The extraction parameters for aperture photometry.

    Returns
    -------
    ra, dec : float
        RA and Dec are the right ascension and declination respectively
        at the nominal center of the image.
    wavelength : ndarray, 1D
        The wavelength in micrometers at each plane of the IFU cube.
    temp_flux : ndarray, 1D
        The sum of the data values in the extraction aperture minus the
        sum of the data values in the background region (scaled by the
        ratio of areas), for each plane.
        The data values are in units of surface brightness, so this value
        isn't really the flux, it's an intermediate value.  Dividing by
        `npixels` (to compute the average) will give the value for the
        `surf_bright` (surface brightness) column, and multiplying by
        the solid angle of a pixel will give the flux for a point source.
    f_var_poisson : ndarray, 1D
        The extracted poisson variance values to go along with the
        temp_flux array.
    f_var_rnoise : ndarray, 1D
        The extracted read noise variance values to go along with the
        temp_flux array.
    f_var_flat : ndarray, 1D
        The extracted flat field variance values to go along with the
        temp_flux array.
    background : ndarray, 1D
        For point source data, the background array is the count rate that was subtracted
        from the total source data values to get `temp_flux`. This background is determined
        for annulus region. For extended source data, the background array is the sigma clipped
        extracted region.
    b_var_poisson : ndarray, 1D
        The extracted poisson variance values to go along with the
        background array.
    b_var_rnoise : ndarray, 1D
        The extracted read noise variance values to go along with the
        background array.
    b_var_flat : ndarray, 1D
        The extracted flat field variance values to go along with the
        background array.
    npixels : ndarray, 1D, float64
        For each slice, this is the number of pixels that were added
        together to get `temp_flux`.
    dq : ndarray, 1D, uint32
        The data quality array.
    npixels_bkg : ndarray, 1D, float64
        For each slice, for point source data  this is the number of pixels that were added
        together to get `temp_flux` for an annulus region or for extended source
        data it is the number of pixels used to determine the background
    radius_match : ndarray,1D, float64
        The size of the extract radius in pixels used at each wavelength of the IFU cube
    x_center, y_center : float
        The x and y center of the extraction region
    """
    data = input_model.data
    try:
        var_poisson = input_model.var_poisson
        var_rnoise = input_model.var_rnoise
        var_flat = input_model.var_flat
    except AttributeError:
        log.info(
            "Input model does not break out variance information. Passing only generalized errors."
        )
        # NIRSpec variances sometimes overflow
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "overflow encountered", RuntimeWarning)
            var_poisson = input_model.err * input_model.err
        var_rnoise = np.zeros_like(input_model.data)
        var_flat = np.zeros_like(input_model.data)
    weightmap = input_model.weightmap

    shape = data.shape
    if len(shape) != 3:
        log.error("Expected a 3D IFU cube; dimension is %d.", len(shape))
        raise RuntimeError("The IFU cube should be 3D.")

    # We need to allocate temp_flux, background, npixels, and dq arrays
    # no matter what.  We may need to divide by npixels, so the default
    # is 1 rather than 0.
    temp_flux = np.zeros(shape[0], dtype=np.float64)
    background = np.zeros(shape[0], dtype=np.float64)
    npixels = np.ones(shape[0], dtype=np.float64)
    npixels_bkg = np.ones(shape[0], dtype=np.float64)
    f_var_poisson = np.zeros(shape[0], dtype=np.float64)
    f_var_rnoise = np.zeros(shape[0], dtype=np.float64)
    f_var_flat = np.zeros(shape[0], dtype=np.float64)
    b_var_poisson = np.zeros(shape[0], dtype=np.float64)
    b_var_rnoise = np.zeros(shape[0], dtype=np.float64)
    b_var_flat = np.zeros(shape[0], dtype=np.float64)

    dq = np.zeros(shape[0], dtype=np.uint32)

    # This boolean mask will be used to mask any bad voxels in a given plane
    bmask = np.zeros([shape[1], shape[2]], dtype=bool)

    x_center, y_center = None, None

    # If the user supplied extraction center coords, use them and
    # ignore all other source type and source position values
    if extract_params["x_center"] is not None:
        x_center = extract_params["x_center"]
        y_center = extract_params["y_center"]
        locn = None

    # For a point source, try to compute the extraction center
    # from the source location
    elif source_type != "EXTENDED":
        # If ifu_autocen is set, try to find the source in the field using DAOphot
        if extract_params["ifu_autocen"] is True:
            log.info("Using auto source detection.")

            # Median collapse across wavelengths, but ignore wavelengths above
            # 26 microns where MRS throughput is extremely low
            (_, _, wavetemp) = get_coordinates(input_model, 1, 1)
            indx = (np.where(wavetemp < 26))[0]
            collapse = np.ma.median(np.ma.masked_invalid(data[indx, :, :]), axis=0)

            # Sigma-clipped stats on collapsed image
            _, clipmed, cliprms = sigclip(collapse)
            # Find source in the collapsed image above 3 sigma
            daofind = DAOStarFinder(fwhm=3.0, threshold=3 * cliprms)
            sources = daofind(collapse - clipmed)
            if sources is None:
                log.warning("Auto source detection failed.")
            else:
                vals = sources["flux"].value
                # Identify brightest source as the target
                indx = np.argmax(vals)
                x_center, y_center = sources[indx]["xcentroid"], sources[indx]["ycentroid"]
                locn = None
                log.info("Auto source detection success.")
                log.info("Using x_center = %g, y_center = %g", x_center, y_center)

        if x_center is None:
            log.info("Using target coordinates.")
            ra_targ = input_model.meta.target.ra
            dec_targ = input_model.meta.target.dec
            locn = locn_from_wcs(input_model, ra_targ, dec_targ)

            if locn is None or np.isnan(locn[0]):
                log.warning(
                    "Couldn't determine source location from WCS, so "
                    "extraction region will be centered."
                )
                x_center = float(shape[-1]) / 2.0
                y_center = float(shape[-2]) / 2.0
            else:
                (x_center, y_center) = locn
                log.info(
                    "Using x_center = %g, y_center = %g, based on TARG_RA and TARG_DEC",
                    x_center,
                    y_center,
                )

    method = extract_params["method"]
    subpixels = extract_params["subpixels"]
    subtract_background = extract_params["subtract_background"]

    width = None
    height = None
    theta = None
    # pull wavelength plane out of input data.
    # using extract 1d wavelength, interpolate the radius,
    # inner_bkg, outer_bkg to match input wavelength

    # find the wavelength array of the IFU cube
    x0 = float(shape[2]) / 2.0
    y0 = float(shape[1]) / 2.0
    (ra, dec, wavelength) = get_coordinates(input_model, x0, y0)

    # interpolate the extraction parameters to the wavelength of the IFU cube
    radius_match = None
    if source_type == "POINT":
        wave_extract = extract_params["wavelength"].flatten()
        inner_bkg = extract_params["inner_bkg"].flatten()
        outer_bkg = extract_params["outer_bkg"].flatten()
        radius = extract_params["radius"].flatten()

        if (input_model.meta.instrument.name == "MIRI") & (
            extract_params["ifu_rscale"] is not None
        ):
            radius = radius * extract_params["ifu_rscale"] / 2.0
            log.info("Scaling radius by factor =  %g", extract_params["ifu_rscale"] / 2.0)

        frad = interp1d(wave_extract, radius, bounds_error=False, fill_value="extrapolate")
        radius_match = frad(wavelength)
        # radius_match is in arc seconds - need to convert to pixels
        # the spatial scale is the same for all wavelengths so we
        # only need to call compute_scale once.

        if locn is None:
            locn_use = (
                input_model.meta.wcsinfo.crval1,
                input_model.meta.wcsinfo.crval2,
                wavelength[0],
            )
        else:
            locn_use = (ra_targ, dec_targ, wavelength[0])

        scale_degrees = compute_scale(
            input_model.meta.wcs, locn_use, disp_axis=input_model.meta.wcsinfo.dispersion_direction
        )

        scale_arcsec = scale_degrees * 3600.00
        radius_match /= scale_arcsec

        finner = interp1d(wave_extract, inner_bkg, bounds_error=False, fill_value="extrapolate")
        inner_bkg_match = finner(wavelength) / scale_arcsec

        fouter = interp1d(wave_extract, outer_bkg, bounds_error=False, fill_value="extrapolate")
        outer_bkg_match = fouter(wavelength) / scale_arcsec

    elif source_type == "EXTENDED":
        # Ignore any input parameters, and extract the whole image.
        width = float(shape[-1])
        height = float(shape[-2])
        x_center = width / 2.0 - 0.5
        y_center = height / 2.0 - 0.5
        theta = 0.0
        subtract_background = False
        bkg_sigma_clip = extract_params["bkg_sigma_clip"]

    log.debug("IFU 1D extraction parameters:")
    log.debug("  x_center = %s", str(x_center))
    log.debug("  y_center = %s", str(y_center))
    if source_type == "POINT":
        log.debug("  method = %s", method)
        if method == "subpixel":
            log.debug("  subpixels = %s", str(subpixels))
    else:
        log.debug("  width = %s", str(width))
        log.debug("  height = %s", str(height))
        log.debug("  theta = %s degrees", str(theta))
        log.debug("  subtract_background = %s", str(subtract_background))
        log.debug("  sigma clip value for background = %s", str(bkg_sigma_clip))
        log.debug("  method = %s", method)
        if method == "subpixel":
            log.debug("  subpixels = %s", str(subpixels))

    position = (x_center, y_center)

    # get aperture for extended it will not change with wavelength
    if source_type == "EXTENDED":
        aperture = RectangularAperture(position, width, height, theta)

    for k in range(shape[0]):  # looping over wavelength
        inner_bkg = None
        outer_bkg = None

        if source_type == "POINT":
            radius = radius_match[k]  # this radius has been converted to pixels
            aperture = CircularAperture(position, r=radius)
            inner_bkg = inner_bkg_match[k]
            outer_bkg = outer_bkg_match[k]
            if inner_bkg <= 0.0 or outer_bkg <= 0.0 or inner_bkg >= outer_bkg:
                log.debug(
                    "Turning background subtraction off, due to "
                    "the values of inner_bkg and outer_bkg."
                )
                subtract_background = False

        if subtract_background and inner_bkg is not None and outer_bkg is not None:
            annulus = CircularAnnulus(position, r_in=inner_bkg, r_out=outer_bkg)
        else:
            annulus = None

        subtract_background_plane = subtract_background
        # Compute the area of the aperture and possibly also of the annulus.
        # for each wavelength bin (taking into account empty spaxels)
        normalization = 1.0
        temp_weightmap = weightmap[k, :, :]
        temp_weightmap[temp_weightmap > 1] = 1
        annulus_area = 0

        # Make a boolean mask to ignore voxels with no valid data
        bmask[:] = False
        bmask[np.where(temp_weightmap == 0)] = True

        # aperture_photometry - using weight map
        phot_table = aperture_photometry(
            temp_weightmap, aperture, mask=bmask, method=method, subpixels=subpixels
        )

        aperture_area = float(phot_table["aperture_sum"][0])
        # if aperture_area = 0, then there is no valid data for this wavelength
        # set the DQ flag to DO_NOT_USE
        if aperture_area == 0:
            dq[k] = dqflags.pixel["DO_NOT_USE"]

        # There is no valid data for this region. To prevent the code from
        # crashing set aperture_area to a nonzero value. It will have the dq flag
        if aperture_area == 0 and aperture.area > 0:
            aperture_area = aperture.area

        if subtract_background and annulus is not None:
            # Compute the area of the annulus.
            phot_table = aperture_photometry(
                temp_weightmap, annulus, mask=bmask, method=method, subpixels=subpixels
            )
            annulus_area = float(phot_table["aperture_sum"][0])

            if annulus_area == 0 and annulus.area > 0:
                annulus_area = annulus.area

            if annulus_area > 0.0:
                normalization = aperture_area / annulus_area
            else:
                log.warning(
                    "Background annulus has no area, so background "
                    f"subtraction will be turned off. {k}"
                )
                subtract_background_plane = False

        npixels[k] = aperture_area

        npixels_bkg[k] = 0.0
        if annulus is not None:
            npixels_bkg[k] = annulus_area
        # aperture_photometry - using data

        phot_table = aperture_photometry(
            data[k, :, :], aperture, mask=bmask, method=method, subpixels=subpixels
        )
        temp_flux[k] = float(phot_table["aperture_sum"][0])

        var_poisson_table = aperture_photometry(
            var_poisson[k, :, :], aperture, mask=bmask, method=method, subpixels=subpixels
        )
        f_var_poisson[k] = float(var_poisson_table["aperture_sum"][0])

        var_rnoise_table = aperture_photometry(
            var_rnoise[k, :, :], aperture, mask=bmask, method=method, subpixels=subpixels
        )
        f_var_rnoise[k] = float(var_rnoise_table["aperture_sum"][0])

        var_flat_table = aperture_photometry(
            var_flat[k, :, :], aperture, mask=bmask, method=method, subpixels=subpixels
        )
        f_var_flat[k] = float(var_flat_table["aperture_sum"][0])

        # Point source type of data with defined annulus size
        if subtract_background_plane:
            bkg_table = aperture_photometry(
                data[k, :, :], annulus, mask=bmask, method=method, subpixels=subpixels
            )
            background[k] = float(bkg_table["aperture_sum"][0])
            temp_flux[k] = temp_flux[k] - background[k] * normalization

            var_poisson_table = aperture_photometry(
                var_poisson[k, :, :], annulus, mask=bmask, method=method, subpixels=subpixels
            )
            b_var_poisson[k] = float(var_poisson_table["aperture_sum"][0])

            var_rnoise_table = aperture_photometry(
                var_rnoise[k, :, :], annulus, mask=bmask, method=method, subpixels=subpixels
            )
            b_var_rnoise[k] = float(var_rnoise_table["aperture_sum"][0])

            var_flat_table = aperture_photometry(
                var_flat[k, :, :], annulus, mask=bmask, method=method, subpixels=subpixels
            )
            b_var_flat[k] = float(var_flat_table["aperture_sum"][0])

            # Propagate errors in background to the background-subtracted science spectrum
            f_var_poisson[k] += b_var_poisson[k] * normalization * normalization
            f_var_rnoise[k] += b_var_rnoise[k] * normalization * normalization
            f_var_flat[k] += b_var_flat[k] * normalization * normalization

        # Extended source data - background determined from sigma clipping
        if source_type == "EXTENDED":
            bkg_data = data[k, :, :]
            # pull out the data with coverage in IFU cube. We do not want to use
            # the edge data that is zero to define the statistics on clipping
            bkg_stat_data = bkg_data[temp_weightmap == 1]

            # If there are good data, work out the statistics
            if len(bkg_stat_data) > 0:
                bkg_mean, _, bkg_stddev = stats.sigma_clipped_stats(
                    bkg_stat_data, sigma=bkg_sigma_clip, maxiters=5
                )
                low = bkg_mean - bkg_sigma_clip * bkg_stddev
                high = bkg_mean + bkg_sigma_clip * bkg_stddev

                # set up the mask to flag data that should not be used in aperture photometry
                # Reject data outside the sigma-clipped range
                maskclip = np.logical_or(bkg_data < low, bkg_data > high)
                # Reject data outside the valid data footprint
                maskclip = np.logical_or(maskclip, bmask)

                bkg_table = aperture_photometry(
                    bkg_data, aperture, mask=maskclip, method=method, subpixels=subpixels
                )
                background[k] = float(bkg_table["aperture_sum"][0])
                phot_table = aperture_photometry(
                    temp_weightmap, aperture, mask=maskclip, method=method, subpixels=subpixels
                )
                npixels_bkg[k] = float(phot_table["aperture_sum"][0])

                var_poisson_table = aperture_photometry(
                    var_poisson[k, :, :],
                    aperture,
                    mask=maskclip,
                    method=method,
                    subpixels=subpixels,
                )
                b_var_poisson[k] = float(var_poisson_table["aperture_sum"][0])

                var_rnoise_table = aperture_photometry(
                    var_rnoise[k, :, :], aperture, mask=maskclip, method=method, subpixels=subpixels
                )
                b_var_rnoise[k] = float(var_rnoise_table["aperture_sum"][0])

                var_flat_table = aperture_photometry(
                    var_flat[k, :, :], aperture, mask=maskclip, method=method, subpixels=subpixels
                )
                b_var_flat[k] = float(var_flat_table["aperture_sum"][0])

        del temp_weightmap
        # done looping over wavelength bins
    # Check for NaNs in the wavelength array, flag them in the dq array,
    # and truncate the arrays if NaNs are found at endpoints (unless the
    # entire array is NaN).

    wavelength, dq, nan_slc = nans_in_wavelength(wavelength, dq)
    temp_flux = temp_flux[nan_slc]
    background = background[nan_slc]
    npixels = npixels[nan_slc]
    npixels_bkg = npixels_bkg[nan_slc]
    f_var_poisson = f_var_poisson[nan_slc]
    f_var_rnoise = f_var_rnoise[nan_slc]
    f_var_flat = f_var_flat[nan_slc]
    b_var_poisson = b_var_poisson[nan_slc]
    b_var_rnoise = b_var_rnoise[nan_slc]
    b_var_flat = b_var_flat[nan_slc]

    return (
        ra,
        dec,
        wavelength,
        temp_flux,
        f_var_poisson,
        f_var_rnoise,
        f_var_flat,
        background,
        b_var_poisson,
        b_var_rnoise,
        b_var_flat,
        npixels,
        dq,
        npixels_bkg,
        radius_match,
        x_center,
        y_center,
    )


def locn_from_wcs(input_model, ra_targ, dec_targ):
    """
    Get the location of the spectrum, based on the WCS.

    Parameters
    ----------
    input_model : JWSTDataModel
        The input science model.
    ra_targ, dec_targ : float or None
        The right ascension and declination of the target (degrees).

    Returns
    -------
    locn : tuple of two int, or None
        X and Y coordinates (in that order) of the pixel that has right
        ascension and declination coordinates closest to the target
        location.  The spectral extraction region should be centered here.
    """
    if ra_targ is None or dec_targ is None:
        log.warning("TARG_RA and/or TARG_DEC not found; can't compute pixel location of target.")
        locn = None
    else:
        shape = input_model.data.shape
        grid = np.indices(shape[-2:])
        z = np.zeros(shape[-2:], dtype=np.float64) + shape[0] // 2
        # The arguments are the X, Y, and Z pixel coordinates, and the
        # output arrays will be 2D.
        (ra_i, dec_i, wl) = input_model.meta.wcs(grid[1], grid[0], z)
        cart = celestial_to_cartesian(ra_i, dec_i)
        v = celestial_to_cartesian(ra_targ, dec_targ)  # a single vector
        diff = cart - v
        # We want the pixel with the minimum distance from v, but the pixel
        # with the minimum value of distance squared will be the same.
        dist2 = (diff**2).sum(axis=-1)
        nan_mask = np.isnan(wl)
        dist2[..., :] = np.where(nan_mask, HUGE_DIST, dist2[..., :])
        del nan_mask
        k = np.argmin(dist2.ravel())
        (j, i) = divmod(k, dist2.shape[1])  # y, x coordinates

        if i <= 0 or j <= 0 or i >= shape[-1] - 1 or j >= shape[-2] - 1:
            log.warning("WCS implies the target is beyond the edge of the image")
            log.warning("This location will not be used")
            locn = None
        else:
            locn = (i, j)  # x, y coordinates

    return locn


def celestial_to_cartesian(ra, dec):
    """
    Convert celestial coordinates to Cartesian.

    Parameters
    ----------
    ra, dec : ndarray or float
        The right ascension and declination (degrees).  Both `ra` and `dec`
        should be arrays of the same shape, or they should both be float.

    Returns
    -------
    cart : ndarray
        If `ra` and `dec` are float, `cart` with be a 3-element array.
        If `ra` and `dec` are arrays, `cart` will be an array with shape
        ra.shape + (3,).
        For each element of `ra` (or `dec`), the last axis of `cart` will
        give the Cartesian coordinates of a unit vector in the direction
        `ra`, `dec`.  The elements of the vector in Cartesian coordinates
        are in the order x, y, z, where x is the direction toward the
        vernal equinox, y is the direction toward right ascension = 90
        degrees (6 hours) and declination = 0, and z is toward the north
        celestial pole.
    """
    if hasattr(ra, "shape"):
        shape = ra.shape + (3,)
    else:
        shape = (3,)

    cart = np.zeros(shape, dtype=np.float64)
    cart[..., 2] = np.sin(dec * np.pi / 180.0)
    r_xy = np.cos(dec * np.pi / 180.0)
    cart[..., 1] = r_xy * np.sin(ra * np.pi / 180.0)
    cart[..., 0] = r_xy * np.cos(ra * np.pi / 180.0)

    return cart


def get_coordinates(input_model, x0, y0):
    """
    Get celestial coordinates and wavelengths.

    Parameters
    ----------
    input_model : IFUCubeModel
        The input model.
    x0, y0 : float
        The pixel number at which the coordinates should be determined.
        If the extract1d reference file is JSON format, this point will be the
        nominal center of the image.  For an image extract1d reference file this
        will be the centroid of the pixels that were flagged as source.

    Returns
    -------
    ra, dec : float
        RA and Dec are the right ascension and declination respectively
        at pixel (0, y0, x0).
    wavelength : ndarray, 1D
        The wavelength in micrometers at each pixel.
    """
    if hasattr(input_model.meta, "wcs"):
        wcs = input_model.meta.wcs
    else:
        log.warning("WCS function not found in input.")
        wcs = None

    nelem = input_model.data.shape[0]

    if wcs is not None:
        x_array = np.empty(nelem, dtype=np.float64)
        x_array.fill(x0)
        y_array = np.empty(nelem, dtype=np.float64)
        y_array.fill(y0)
        z_array = np.arange(nelem, dtype=np.float64)  # for wavelengths
        ra, dec, wavelength = wcs(x_array, y_array, z_array)
        # ra and dec should be constant, so nelem // 2 is an arbitrary element.
        ra = ra[nelem // 2]
        dec = dec[nelem // 2]
    else:
        (ra, dec) = (0.0, 0.0)
        wavelength = np.arange(1, nelem + 1, dtype=np.float64)

    return ra, dec, wavelength


def nans_in_wavelength(wavelength, dq):
    """
    Check for NaNs in the wavelength array.

    If NaNs are found in the wavelength array, flag them in the dq array,
    and truncate the arrays at either or both ends if NaNs are found at
    endpoints (unless the entire array is NaN).

    Parameters
    ----------
    wavelength : ndarray, 1D, float64
        The wavelength in micrometers at each pixel.
    dq : ndarray, 1D, uint32
        The data quality array.

    Returns
    -------
    wavelength : ndarray, 1D, float64
        Truncated wavelength array
    dq : ndarray, 1D, uint32
        Truncated DQ array.
    slc : slice
        Slice used for truncation.
    """
    nelem = np.size(wavelength)
    slc = slice(nelem)
    if nelem == 0:
        log.warning("Output arrays are empty!")
        return wavelength, dq, slice(nelem)

    nan_mask = np.isnan(wavelength)
    n_nan = nan_mask.sum(dtype=np.intp)
    if n_nan == nelem:
        log.warning("Wavelength array is all NaN!")
        dq = np.bitwise_or(dq[:], dqflags.pixel["DO_NOT_USE"])
        return wavelength, dq, slice(0)

    if n_nan > 0:
        log.warning("%d NaNs in wavelength array.", n_nan)
        dq[nan_mask] = np.bitwise_or(dq[nan_mask], dqflags.pixel["DO_NOT_USE"])
        not_nan = np.logical_not(nan_mask)
        flag = np.where(not_nan)
        if len(flag[0]) > 0:
            n_trimmed = flag[0][0] + nelem - (flag[0][-1] + 1)
            if n_trimmed > 0:
                slc = slice(flag[0][0], flag[0][-1] + 1)
                wavelength = wavelength[slc]
                dq = dq[slc]
                log.info("Output arrays have been trimmed by %d elements", n_trimmed)

    return wavelength, dq, slc


def separate_target_and_background(ref):
    """
    Create masks for target and background.

    Parameters
    ----------
    ref : ndarray, 2D or 3D
        This is the reference image data array.  This should be the same
        shape as one plane of the science data (or the same shape as the
        entire 3D science data array).  A value of 1 in a pixel indicates
        that the pixel should be included when computing the target
        spectrum.  A value of 0 means the pixel is not part of either the
        target or background.  A value of -1 means the pixel should be
        included as part of the background region.

    Returns
    -------
    mask_target : ndarray, 2D or 3D
        This is an array of the same type and shape as the reference
        image, but with values of only 0 or 1.  A value of 1 indicates
        that the corresponding pixel of the science data array should be
        included when adding up values to make the 1D spectrum, and a
        value of 0 means that it should not be included.
    mask_bkg : ndarray, 2D or 3D, or None.
        This is like `mask_target` (i.e. values are 0 or 1) but for
        background regions.  A value of -1 in `mask_bkg` indicates a pixel
        that should be included as part of the background.  If there is no
        pixel in the reference image with a value of -1, `mask_bkg` will
        be set to None.
    """
    mask_target = np.where(ref == 1.0, 1.0, 0.0)

    if np.any(ref == -1.0):
        mask_bkg = np.where(ref == -1.0, 1.0, 0.0)
    else:
        mask_bkg = None

    return mask_target, mask_bkg


def im_centroid(data, mask_target):
    """
    Compute the mean location of the target.

    Parameters
    ----------
    data : ndarray, 3D
        This is the science image data array.
    mask_target : ndarray, 2D or 3D
        This is an array of the same type and shape as one plane of the
        science image (or the same type and shape of the entire 3D science
        image), but with values of 0 or 1, where 1 indicates a pixel within
        the source.

    Returns
    -------
    y0, x0 : tuple of two float
        The centroid of pixels flagged as source.
    """
    # Collapse the science data along the dispersion direction to get a
    # 2D image of the IFU field of view.  Multiplying by mask_target
    # zeros out all pixels that are not regarded as part of the target
    # (or targets).
    if len(mask_target.shape) == 2:
        data_2d = data.sum(axis=0, dtype=np.float64) * mask_target
    else:
        data_2d = (data * mask_target).sum(axis=0, dtype=np.float64)
    if data_2d.sum() == 0.0:
        log.warning("Couldn't compute image centroid.")
        shape = data_2d.shape
        y0 = shape[0] / 2.0
        x0 = shape[0] / 2.0
        return y0, x0

    x_profile = data_2d.sum(axis=0, dtype=np.float64)
    x = np.arange(data_2d.shape[1], dtype=np.float64)
    s_x = (x_profile * x).sum()
    x0 = s_x / x_profile.sum()

    y_profile = data_2d.sum(axis=1, dtype=np.float64)
    y = np.arange(data_2d.shape[0], dtype=np.float64)
    s_y = (y_profile * y).sum()
    y0 = s_y / y_profile.sum()

    return y0, x0


def shift_ref_image(mask, delta_y, delta_x, fill=0):
    """
    Apply source position offset to target or background for ref image.

    Parameters
    ----------
    mask : ndarray, 2D or 3D
        This is either the target mask or the background mask, which was
        created from an image extract1d reference file.
        It is assumed that all pixels in the science data are good.  If
        that is not correct, `shift_ref_image` may be called with a data
        quality array as the first argument, in order to obtain a shifted
        data quality array.  In this case, `fill` should be set to a
        positive value, e.g. 1.
    delta_y, delta_x : int
        These are the shifts to be applied to the vertical and horizontal
        axes respectively.
    fill : int or float
        The output array will be initialized to this value.  This should
        be 0 (the default) if `mask` is a target or background mask, but
        it should be set to 1 (or some other positive value) if `mask` is
        a data quality array.

    Returns
    -------
    temp : ndarray, same type and shape as `mask`
        A copy of `mask`, but shifted by `delta_y` and `delta_x`.
    """
    if delta_x == 0 and delta_y == 0:
        return mask.copy()

    shape = mask.shape
    if abs(delta_y) >= shape[-2] or abs(delta_x) >= shape[-1]:
        log.warning("Nod offset %d or %d is too large, skipping ...", delta_y, delta_x)
        return None

    if delta_y > 0:
        islice_y = slice(0, -delta_y)
        oslice_y = slice(delta_y, shape[-2])
    elif delta_y < 0:
        islice_y = slice(-delta_y, shape[-2])
        oslice_y = slice(0, delta_y)
    else:
        islice_y = slice(0, shape[-2])
        oslice_y = slice(0, shape[-2])

    if delta_x > 0:
        islice_x = slice(0, -delta_x)
        oslice_x = slice(delta_x, shape[-1])
    elif delta_x < 0:
        islice_x = slice(-delta_x, shape[-1])
        oslice_x = slice(0, delta_x)
    else:
        islice_x = slice(0, shape[-1])
        oslice_x = slice(0, shape[-1])

    temp = np.zeros_like(mask) + fill
    temp[..., oslice_y, oslice_x] = mask[..., islice_y, islice_x]

    return temp


def sigma_clip_extended_region(
    data, var_poisson, var_rnoise, var_flat, mask_targ, wmap, sigma_clip
):
    """
    Sigma clip the extraction region.

    Parameters
    ----------
    data : ndarray, 3D
        Input data array to perform extraction from.
    var_poisson : ndarray, 2D
        Poisson noise variance array to be extracted following data extraction method.
    var_rnoise : ndarray, 2D
        Read noise variance array to be extracted following data extraction method.
    var_flat : ndarray, 2D
        Flat noise variance array to be extracted following data extraction method.
    mask_targ : ndarray, 2D or 3D
        Mask of pixels defining the extended source region. A value of 1 indicated
        pixel is in the extraction region.
    wmap : ndarray, 3D
        Weight map for IFU.
    sigma_clip : float
        Outlier sigma clipping parameter.

    Returns
    -------
    sigma_clip_region : ndarray, 1D
        Summed extracted region with sigma clipping for each wavelength plane.
    d_var_poisson : ndarray, 1D
        Sigma-clipped var_poisson array.
    d_var_rnoise : ndarray, 1D
        Sigma-clipped var_rnoise array.
    d_var_flat : ndarray, 1D
        Sigma-clipped var_flat array.
    n_bkg : ndarray, 1D
        Sum of pixels used in sigma clipped extracted region.
    """
    shape = data.shape
    shape_ref = mask_targ.shape
    n_bkg = np.ones(shape[0], dtype=np.float64)

    sigma_clip_region = np.zeros(shape[0], dtype=np.float64)
    d_var_poisson = np.zeros(shape[0], dtype=np.float64)
    d_var_rnoise = np.zeros(shape[0], dtype=np.float64)
    d_var_flat = np.zeros(shape[0], dtype=np.float64)

    # for each wavelength plane mark outliers as 0 in mask_bkg
    for k in range(shape[0]):  # looping over wavelength
        if len(shape_ref) == 2:
            extract_region = mask_targ.copy()
        else:
            extract_region = mask_targ[k, :, :].copy()
        data_plane = data[k, :, :]
        var_poisson_plane = var_poisson[k, :, :]
        var_rnoise_plane = var_rnoise[k, :, :]
        var_flat_plane = var_flat[k, :, :]
        # pull out extract source region to determined stats on for sigma clipping
        extract_data = data_plane[extract_region == 1]
        ext_mean, _, ext_stddev = stats.sigma_clipped_stats(
            extract_data, sigma=sigma_clip, maxiters=5
        )
        low = ext_mean - sigma_clip * ext_stddev
        high = ext_mean + sigma_clip * ext_stddev

        # set up the mask to flag data that should not be used
        maskclip = np.logical_or(data_plane < low, data_plane > high)  # flag outliers
        maskclip = np.logical_or(maskclip, ~np.isfinite(data_plane))
        extract_region[maskclip] = 0

        sigma_clip_region[k] = np.sum(data_plane * extract_region * wmap[k, :, :])
        d_var_poisson[k] = np.sum(var_poisson_plane * extract_region * wmap[k, :, :])
        d_var_rnoise[k] = np.sum(var_rnoise_plane * extract_region * wmap[k, :, :])
        d_var_flat[k] = np.sum(var_flat_plane * extract_region * wmap[k, :, :])
        n_bkg[k] = np.sum(wmap[k, :, :] * extract_region)

    return sigma_clip_region, d_var_poisson, d_var_rnoise, d_var_flat, n_bkg
