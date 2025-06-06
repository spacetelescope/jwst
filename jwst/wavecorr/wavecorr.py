import logging
import numpy as np
from gwcs import wcstools
from gwcs import coordinate_frames as cf
from astropy import units as u
from astropy.modeling import tabular
from astropy.modeling.mappings import Identity

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, wavecorr_file):
    """
    Perform wavelength correction for NIRSpec MOS and FS point sources.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    wavecorr_file : str
        Wavecorr reference file name.

    Returns
    -------
    output_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Corrected data model.
    """
    wavecorr_supported_modes = ["NRS_FIXEDSLIT", "NRS_MSASPEC", "NRS_BRIGHTOBJ", "NRS_AUTOFLAT"]

    # Check for valid exposure type
    exp_type = input_model.meta.exposure.type.upper()
    if exp_type not in wavecorr_supported_modes:
        log.info(f"Skipping wavecorr correction for EXP_TYPE {exp_type}")
        input_model.meta.cal_step.wavecorr = "SKIPPED"
        return input_model

    output_model = input_model.copy()

    # For BRIGHTOBJ, operate on the single SlitModel
    corrected = False
    if isinstance(input_model, datamodels.SlitModel):
        if _is_point_source(input_model):
            corrected = apply_zero_point_correction(output_model, wavecorr_file)
    else:
        for slit in output_model.slits:
            if _is_point_source(slit):
                completed = apply_zero_point_correction(slit, wavecorr_file)
                if completed:
                    corrected = True
                    slit.wavelength_corrected = True
                else:  # pragma: no cover
                    log.warning(f"Corrections are not invertible for slit {slit.name}")
                    log.warning("Skipping wavecorr correction")
                    slit.wavelength_corrected = False
            else:
                slit.wavelength_corrected = False

    if corrected:
        output_model.meta.cal_step.wavecorr = "COMPLETE"
    else:
        output_model.meta.cal_step.wavecorr = "SKIPPED"

    return output_model


def apply_zero_point_correction(slit, reffile):
    """
    Apply the NIRSpec wavelength zero-point correction.

    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`, `~jwst.datamodels.CubeModel`
        Slit data to be corrected.
    reffile : str
        The ``wavecorr`` reference file.

    Returns
    -------
    completed : bool
        A flag to report whether the zero-point correction was added or skipped.
    """
    slit_wcs = slit.meta.wcs

    # Retrieve the source position and aperture name from metadata
    source_xpos = slit.source_xpos
    if slit.meta.exposure.type in ["NRS_FIXEDSLIT", "NRS_BRIGHTOBJ"]:
        aperture_name = slit.name
    else:
        # For the MSA the aperture name is "MOS"
        aperture_name = "MOS"

    lam = slit.wavelength.copy() * 1e-6
    dispersion = compute_dispersion(slit.meta.wcs)

    wave2wavecorr = calculate_wavelength_correction_transform(
        lam, dispersion, reffile, source_xpos, aperture_name
    )

    # wave2wavecorr should not be None for real data
    if wave2wavecorr is None:
        completed = False
        return completed
    else:
        # Make a new frame to insert into the slit wcs object
        slit_spatial = cf.Frame2D(
            name="slit_spatial", axes_order=(0, 1), unit=("", ""), axes_names=("x_slit", "y_slit")
        )
        spec = cf.SpectralFrame(
            name="spectral", axes_order=(2,), unit=(u.micron,), axes_names=("wavelength",)
        )
        wcorr_frame = cf.CompositeFrame([slit_spatial, spec], name="wavecorr_frame")

        # Insert the new transform into the slit wcs object
        wave2wavecorr = Identity(2) & wave2wavecorr
        slit_wcs.insert_frame("slit_frame", wave2wavecorr, wcorr_frame)

        # Update the stored wavelengths for the slit
        slit.wavelength = compute_wavelength(slit_wcs)

        completed = True
        return completed


def calculate_wavelength_correction_transform(
    lam, dispersion, freference, source_xpos, aperture_name
):
    """
    Add a transform for wavelength correction to the WCS pipeline for each slit.

    Parameters
    ----------
    lam : ndarray
        Wavelength array [in m].
    dispersion : ndarray
        The pixel dispersion [in m].
    freference : str
        ``wavecorr`` reference file name.
    source_xpos : float
        X position of the source as a fraction of the slit size.
    aperture_name : str
        Aperture name.

    Returns
    -------
    model : `~astropy.modeling.tabular.Tabular1D`or None
        A model which takes wavelength inputs and returns zero-point
        corrected wavelengths.  Returns None if an invertible model
        cannot be generated.
    """
    # Open the zero point reference model
    with datamodels.WaveCorrModel(freference) as wavecorr:
        for ap in wavecorr.apertures:
            if ap.aperture_name == aperture_name:
                log.info(f"Using wavelength zero-point correction for aperture {ap.aperture_name}")
                offset_model = ap.zero_point_offset.copy()
                break
        else:
            log.warning(f"No wavelength zero-point correction found for slit {aperture_name}")
            return None

    # Set lookup table to extrapolate at bounds to recover wavelengths
    # beyond model bounds, particularly for the red and blue ends of
    # prism observations.  fill_value = None sets the lookup tables
    # to use the default extrapolation which is a linear extrapolation
    # from scipy.interpolate.interpn
    offset_model.bounds_error = False
    offset_model.fill_value = None

    # Average the wavelength and dispersion across 2D extracted slit and remove nans
    # So that we have a 1D wavelength array for building a 1D lookup table wcs transform
    lam_mean = np.nanmean(lam, axis=0)
    disp_mean = np.nanmean(dispersion, axis=0)
    nan_lams = np.isnan(lam_mean) | np.isnan(disp_mean)
    lam_mean = lam_mean[~nan_lams]
    disp_mean = disp_mean[~nan_lams]

    # Calculate the corrected wavelengths
    pixel_corrections = offset_model(lam_mean, source_xpos)
    lam_corrected = lam_mean + (pixel_corrections * disp_mean)

    # Check to make sure that the corrected wavelengths are monotonically increasing
    if np.all(np.diff(lam_corrected) > 0):
        # monotonically increasing
        # Build a look up table to transform between corrected and uncorrected wavelengths
        wave2wavecorr = tabular.Tabular1D(
            points=lam_mean,
            lookup_table=lam_corrected,
            bounds_error=False,
            fill_value=None,
            name="wave2wavecorr",
        )

        return wave2wavecorr

    else:
        # output wavelengths are not monotonically increasing
        return None


def compute_dispersion(wcs, xpix=None, ypix=None):
    """
    Compute the pixel dispersion.

    If `xpix` or `ypix` is not provided, the dispersion is computed on a grid
    based on ``wcs.bounding_box``.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    xpix : ndarray, float, optional
        X pixel coordinates for dispersion elements.
    ypix : ndarray, float, optional
        Y pixel coordinates for dispersion elements.

    Returns
    -------
    dispersion : ndarray
        The pixel dispersion [in m].
    """
    if xpix is None or ypix is None:
        xpix, ypix = wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1, 1))
    xleft = xpix - 0.5
    xright = xpix + 0.5
    _, _, lamright = wcs(xright, ypix)
    _, _, lamleft = wcs(xleft, ypix)
    return (lamright - lamleft) * 10**-6


def compute_wavelength(wcs, xpix=None, ypix=None):
    """
    Compute the pixel wavelength.

    If `xpix` or `ypix` is not provided, the dispersion is computed on a grid
    based on ``wcs.bounding_box``.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    xpix : ndarray, float, optional
        X pixel coordinates for dispersion elements.
    ypix : ndarray, float, optional
        Y pixel coordinates for dispersion elements.

    Returns
    -------
    wavelength : ndarray
        The wavelength [in microns].
    """
    if xpix is None or ypix is None:
        xpix, ypix = wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1, 1))

    _, _, lam = wcs(xpix, ypix)
    return lam


def _is_point_source(slit):
    """
    Determine if a source is a point source.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.transforms.models.Slit`
        A slit object.

    Returns
    -------
    bool
        True if point source; False otherwise.
    """
    result = False

    # Get the source type value set by the source_type step (if any)
    if slit.source_type is not None:
        src_type = slit.source_type
    elif slit.meta.target.source_type is not None:
        src_type = slit.meta.target.source_type
    else:
        src_type = None

    if src_type is not None and src_type.upper() in ["POINT", "EXTENDED"]:
        # Use the supplied value
        log.info(f"Detected a {src_type} source type in slit {slit.name}")
        if src_type.strip().upper() == "POINT":
            result = True
        else:
            result = False
    else:
        log.info("Unknown source type")

    return result
