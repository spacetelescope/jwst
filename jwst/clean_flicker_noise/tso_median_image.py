import logging
import warnings

import numpy as np
from stdatamodels.jwst import datamodels

from jwst.background.background_step import BackgroundStep
from jwst.datamodels.utils.tso_multispec import make_tso_specmodel
from jwst.extract_1d.extract_1d_step import Extract1dStep
from jwst.extract_1d.soss_extract import pastasoss, soss_extract
from jwst.extract_2d.extract_2d_step import Extract2dStep
from jwst.lib.basic_utils import disable_logging
from jwst.tso_photometry.tso_photometry_step import TSOPhotometryStep
from jwst.white_light.white_light_step import WhiteLightStep

__all__ = ["make_median_image"]

log = logging.getLogger(__name__)


def _soss_box_extract(rateints, soss_refmodel=None):
    """
    Extract spectra with a simple box around the SOSS trace.

    The spectra are intended for approximate scaling, so they do
    not need to be very accurate.  Only order 0 is extracted from
    each integration, with a width of 15 pixels around the trace from
    the pastasoss model.

    Parameters
    ----------
    rateints : `~stdatamodels.jwst.datamodels.CubeModel`
        Background subtracted rateints datamodel.
    soss_refmodel : `~stdatamodels.jwst.datamodels.PastasossModel`, optional
        Used to identify SOSS traces for box extraction.
        If not provided, a default model will be retrieved.

    Returns
    -------
    multi_spec : ~stdatamodels.jwst.datamodels.TSOMultiSpecModel`
        Extracted spectra, with only FLUX and WAVELENGTH arrays
        populated.
    """
    nints = rateints.data.shape[0]
    img_shape = rateints.data.shape[-2:]
    pwcpos = rateints.meta.instrument.pupil_position
    subarray = rateints.meta.subarray.name

    if soss_refmodel is None:
        soss_refmodel = pastasoss.retrieve_default_pastasoss_model()

    # Set a reasonable extraction width
    width = 15

    # Extract only order 1 for scaling purposes
    order = 1

    # Make a multispec model to hold the spectra
    spec_list = []
    for i in range(nints):
        sci_data = rateints.data[i]
        sci_mask = (rateints.dq[i] & datamodels.dqflags.pixel["DO_NOT_USE"]) > 0

        # Get trace x,y positions
        _, xtrace, ytrace, _ = pastasoss.get_soss_traces(
            pwcpos, order=order, refmodel=soss_refmodel
        )
        box_weights = soss_extract.get_box_weights(
            ytrace, width, img_shape, cols=xtrace.astype(int)
        )

        wavemaps = pastasoss.get_soss_wavemaps(
            pwcpos, subarray=subarray, refmodel=soss_refmodel, orders_requested=[1]
        )
        wave_grid = wavemaps[0]

        box_weights[sci_mask] = 0
        npix = np.sum(box_weights, axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            flux = np.nansum(sci_data * box_weights, axis=0)
            wavelength = np.nansum(wave_grid * box_weights, axis=0) / npix

        # Store fluxes with valid wavelengths
        valid = ~np.isnan(wavelength)
        spec = datamodels.SpecModel()
        spec.spec_table = np.zeros((valid.sum(),), dtype=spec.spec_table.dtype)
        spec.spec_table["FLUX"] = flux[valid]
        spec.spec_table["WAVELENGTH"] = wavelength[valid]
        spec_list.append(spec)

    tso_spec = make_tso_specmodel(spec_list)

    # populate the midtime array with unique placeholder values to
    # make sure every integration appears in the whitelight table, later
    tso_spec.spec_table["MJD-AVG"] = np.arange(nints)

    multi_spec = datamodels.TSOMultiSpecModel()
    multi_spec.update(rateints)
    multi_spec.spec.append(tso_spec)
    return multi_spec


def _make_background_ramp(input_ramp, background_rateints):
    """
    Extrapolate a background rate to an up-the-ramp sampling.

    Given a ramp and a model of the background rate per integration,
    produce an up-the-ramp version of the countrate image to match the ramp.

    Total time for one group is (nframes + groupgap) * frame_time. The group
    gap occurs after the last frame, so the read time for group i
    (zero-indexed) is:

    group_time = (nframes + (nframes + groupgap) * i) * frame_time

    Parameters
    ----------
    input_ramp : `~stdatamodels.jwst.datamodels.RampModel`
        The input ramp datamodel to match.
    background_rateints : ndarray
        A 3D background array, with one rate image per integration.

    Returns
    -------
    background_ramp : ndarray
        4D up-the-ramp background data, matching the array shape of
        the input ramp.
    """
    nints, ngroups, ny, nx = input_ramp.data.shape
    nframes = input_ramp.meta.exposure.nframes
    groupgap = input_ramp.meta.exposure.groupgap
    frame_time = input_ramp.meta.exposure.frame_time
    one_group = nframes + groupgap

    background_ramp = np.zeros_like(input_ramp.data)
    for i in range(ngroups):
        group_time = (nframes + one_group * i) * frame_time
        background_ramp[:, i, ...] = background_rateints * group_time

    return background_ramp


def make_median_image(input_model, rateints_model, soss_refmodel=None):
    """
    Make a scaled median image across integrations to subtract before cleaning.

    The procedure is:

    1. Subtract a reference background from the input ramp and rate data
       (NIRISS SOSS only).
    2. Compute a representative flux for scaling for each integration
       from the rate data.
       a. For spectral modes, extract a spectrum from a simple box
          and compute the whitelight flux for each integration.
       b. For imaging modes, compute the overall median of the image.
    3. Median combine all data across integrations to make a median
       ramp or image.
    4. Scale the median image by the representative flux for each
       integration.
    5. Add the subtracted background level to the median image
       (NIRISS SOSS only).
    6. Replace any NaN values in the median image with zeros.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.RampModel` or \
                  `~stdatamodels.jwst.datamodels.CubeModel`
        Ramp or rateints model to be cleaned.
    rateints_model : `~stdatamodels.jwst.datamodels.CubeModel`
        Draft rateints model corresponding to the input model.
        May be the same as ``input_model``.
    soss_refmodel : `~stdatamodels.jwst.datamodels.PastasossModel`, optional
        Used to identify NIRISS SOSS traces for box extraction.
        If not provided, a default model will be retrieved.

    Returns
    -------
    median_image : ndarray
        The scaled median image to subtract, matching the dimensions of
        the input model.

    Raises
    ------
    ValueError
        If the input does not have multiple integrations.
    """
    ndim = input_model.data.ndim
    if ndim < 3:
        raise ValueError("Cannot make a median image for 2D data")
    nints = input_model.data.shape[0]
    if nints <= 1:
        raise ValueError("Cannot make a median image for <2 integrations")

    # Run background subtraction for rateints files
    exp_type = input_model.meta.exposure.type
    if exp_type == "NIS_SOSS":
        log.info("Subtracting SOSS background")
        with disable_logging(level=logging.WARNING):
            step = BackgroundStep()
            bgsub_rateints = step.run(rateints_model)

        # Subtract rateints models to get the background by integration
        background_rate = rateints_model.data - bgsub_rateints.data

    else:
        bgsub_rateints = rateints_model
        background_rate = 0.0

    # Box extraction for flux scaling
    if exp_type == "NIS_SOSS":
        # Simple direct box extraction for SOSS
        multi_spec = _soss_box_extract(bgsub_rateints, soss_refmodel=soss_refmodel)
    elif exp_type in ["NRS_BRIGHTOBJ", "NRC_TSGRISM"]:
        with disable_logging(level=logging.WARNING):
            # Run extract2d to assign a slit-appropriate WCS
            # (required for extract_1d)
            step = Extract2dStep()
            single_slit = step.run(bgsub_rateints)

            # Set the source type to POINT for TSO
            single_slit.source_type = "POINT"

            # Call extract_1d with latest CRDS parameters (via call)
            # since the recommended defaults may vary by mode
            multi_spec = Extract1dStep.call(single_slit, save_results=False)
            single_slit.close()
    else:
        # Not a TSO spectral mode
        multi_spec = None

    # Sum the flux for each integration.
    # Use the whitelight step to retrieve and use a wavelengthrange
    # file as appropriate.
    if multi_spec is not None:
        log.info("Computing an approximate whitelight curve for scaling")
        with disable_logging(level=logging.WARNING):
            step = WhiteLightStep()
            whitelight_table = step.run(multi_spec)
        multi_spec.close()

        if exp_type == "NRS_BRIGHTOBJ":
            detector = input_model.meta.instrument.detector
            wlc_flux = whitelight_table[f"whitelight_flux_{detector}"]
        else:
            wlc_flux = whitelight_table["whitelight_flux"]
        norm_flux = wlc_flux / np.nanmedian(wlc_flux)

    else:
        # Call tso_photometry with latest CRDS parameters (via call),
        # since the recommended defaults may vary by mode.
        with disable_logging(level=logging.WARNING):
            phot_table = TSOPhotometryStep.call(bgsub_rateints)

        # Use the aperture sum as the representative flux for scaling
        phot_flux = phot_table["aperture_sum"].value
        norm_flux = phot_flux / np.nanmedian(phot_flux)

    # Make a background corrected ramp
    if ndim > 3:
        background_ramp = _make_background_ramp(input_model, background_rate)
        bgsub_ramp = input_model.data - background_ramp
    else:
        background_ramp = background_rate
        bgsub_ramp = bgsub_rateints.data

    # Make a median background-subtracted ramp
    log.info("Making a scaled median image")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median_ramp = np.nanmedian(bgsub_ramp, axis=0)

    # Scale the median ramp by the normalized flux
    if ndim == 3:
        scaled_median_ramp = norm_flux[:, None, None] * median_ramp[None, ...]
    else:
        scaled_median_ramp = norm_flux[:, None, None, None] * median_ramp[None, ...]

    # Final data to subtract is background ramp + scaled median
    median_image = background_ramp + scaled_median_ramp

    # Set any NaN pixels to zero for subtraction
    median_image[~np.isfinite(median_image)] = 0.0

    return median_image
