#
#  Module for applying wavelength corrections to NIRSpec MOS and FS
#  slit instances in which the source is off center.
#
import logging
import numpy as np
from gwcs import wcstools

from .. import datamodels
from ..transforms import models as trmodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, wavecorr_file):
    """
    Main wavecorr correction function for NIRSpec MOS and FS.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    wavecorr_file : str
        Wavecorr reference file name.
    """

    wavecorr_supported_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                                'NRS_AUTOFLAT']

    # Check for valid exposure type
    exp_type = input_model.meta.exposure.type.upper()
    if exp_type not in wavecorr_supported_modes:
        log.info(f'Skipping wavecorr correction for EXP_TYPE {exp_type}')
        return input_model

    output_model = input_model.copy()

    # Get the primary slit for a FS exposure
    if exp_type == 'NRS_FIXEDSLIT':
        primary_slit = input_model.meta.instrument.fixed_slit
        if primary_slit is None or primary_slit == 'NONE':
            log.warning('Primary slit name not found in input')
            log.warning('Skipping wavecorr correction')
            input_model.meta.cal_step.wavecorr = 'SKIPPED'
            return input_model

    # For BRIGHTOBJ, operate on the single SlitModel
    if isinstance(input_model, datamodels.SlitModel):
        if _is_point_source(input_model, exp_type):
            apply_zero_point_correction(output_model, wavecorr_file)

    else:

        # For FS only work on the primary slit
        if exp_type == 'NRS_FIXEDSLIT':
            for slit in output_model.slits:
                if slit.name == primary_slit and _is_point_source(slit, exp_type):
                    apply_zero_point_correction(slit, wavecorr_file)
                    break

        # For MOS work on all slits containing a point source
        else:
            for slit in output_model.slits:
                if _is_point_source(slit, exp_type):
                    apply_zero_point_correction(slit, wavecorr_file)

    output_model.meta.cal_step.wavecorr = 'COMPLETE'

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
    """
    slit_wcs = slit.meta.wcs

    # Get the source position in the slit and set the aperture name
    if slit.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        # pass lam = 2 microns
        # needed for wavecorr with fixed slits
        msa_model = get_msa_model(slit)
        source_xpos = get_source_xpos(slit, slit_wcs, lam=2,
                                      msa_model=msa_model)
        aperture_name = slit.name
    else:
        source_xpos = slit.source_xpos
        # For the MSA the aperture name is "MOS"
        aperture_name = "MOS"

    lam = slit.wavelength
    dispersion = compute_dispersion(slit.meta.wcs)
    corr, dq_lam = compute_zero_point_correction(lam, reffile, source_xpos,
                                                 aperture_name, dispersion)
    ## TODO: set a DQ flag to a TBD value for pixels where dq_lam == 0.
    ## The only purpose of dq_lam is to set that flag.

    # Wavelength is in um, the correction is computed in meters.
    slit.wavelength = lam - corr * 10 ** 6


def compute_zero_point_correction(lam, freference, source_xpos, aperture_name, dispersion):
    """
    Compute the NIRSpec wavelength zero-point correction.

    Parameters
    ----------
    lam : ndarray
        Wavelength array.
    freference : str
        ``wavecorr`` reference file name.
    source_xpos : float
        X position of the source as a fraction of the slit size.
    aperture_name : str
        Aperture name.
    dispersion : ndarray
        The pixel dispersion [in m].

    Returns
    -------
    lambda_corr : ndarray
        Wavelength correction.
    lam : ndarray
        Interpolated wavelengths. Extrapolated values are reset to 0.
        This is returned so that the DQ array can be updated with a flag
        which indicates that no zero-point correction was done.
    """
    with datamodels.WaveCorrModel(freference) as wavecorr:
        for ap in wavecorr.apertures:
            if ap.aperture_name == aperture_name:
                log.info(f'Using wavelength zero-point correction for aperture {ap.aperture_name}')
                offset_model = ap.zero_point_offset.copy()
                # TODO: implement variance
                #variance = ap.variance.copy()
                # width = ap.width
                break
        else:
            log.info(f'No wavelength zero-point correction found for slit {aperture_name}')

    deltax = source_xpos
    lam = lam.copy()
    l = lam[~np.isnan(lam)]
    offset_model.bounds_error = False
    correction = offset_model(l * 10 ** -6, [deltax] * l.size)
    lam[~np.isnan(lam)] = correction

    # The correction for pixels outside the slit and wavelengths
    # outside the wave_range is 0.
    lam[np.isnan(lam)] = 0.
    lambda_cor = dispersion * lam
    return lambda_cor, lam


def compute_dispersion(wcs, xpix=None, ypix=None):
    """
    Compute the pixel dispersion.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    xpix : ndarray, float, optional
    ypix : ndarray, float, optional
        Compute the dispersion at the x, y pixels.
        If not provided the dispersion is computed on a
        grid based on ``wcs.bounding_box``.

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
    return (lamright - lamleft) * 10 ** -6


def _is_point_source(slit, exp_type):
    """
    Determine if a source is a point source.

    Parameters
    ----------
    slit : `~jwst.transforms.models.Slit`
        A slit object.
    exp_type : str
        The exposure type
    """
    result = False

    # Get the source type value set by the source_type step (if any)
    if slit.source_type is not None:
        src_type = slit.source_type
    elif slit.meta.target.source_type is not None:
        src_type = slit.meta.target.source_type
    else:
        src_type = None

    if (src_type is not None) and (src_type.upper() in ['POINT', 'EXTENDED']):
        # Use the supplied value
        log.info(f'Detected a {src_type} source type in slit {slit.name}')
        if src_type.strip().upper() == 'POINT':
            result = True
        else:
            result = False
    else:
        log.info("Unknown source type")

    return result


def get_source_xpos(slit, slit_wcs, lam, msa_model):
    """
    Compute the source position within the slit for a NIRSpec fixed slit.

    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`
        The slit object.
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    lam : float
        Wavelength in microns.
    msa_model : `~jwst.datamodels.MSAModel`
        NIRSpec MSA model.

    Returns
    -------
    xpos : float
        X coordinate of the source as a fraction of the slit size.
    """
    xoffset = slit.meta.dither.x_offset  # in arcsec
    yoffset = slit.meta.dither.y_offset  # in arcsec
    v2ref = slit.meta.wcsinfo.v2_ref  # in arcsec
    v3ref = slit.meta.wcsinfo.v3_ref  # in arcsec
    v3idlyangle = slit.meta.wcsinfo.v3yangle  # in deg
    vparity = slit.meta.wcsinfo.vparity

    idl2v23 = trmodels.IdealToV2V3(v3idlyangle, v2ref, v3ref, vparity)
    # Compute the location in V2,V3 [in arcsec]
    xv, yv = idl2v23(xoffset, yoffset)

    v2v3_to_msa_frame = slit_wcs.get_transform("v2v3", "msa_frame")
    xpos_abs, ypos_abs, lam = v2v3_to_msa_frame(xv, yv, lam)
    xpos_frac = absolute2fractional(msa_model, slit, xpos_abs, ypos_abs)
    return xpos_frac


def get_msa_model(input_model):
    """Get the reference file used in constructing the WCS.
    """
    msa_ref = input_model.meta.ref_file.msa.name
    from .. import assign_wcs
    from .. datamodels import MSAModel
    step = assign_wcs.AssignWcsStep()
    msa = MSAModel(step.reference_uri_to_cache_path(msa_ref))
    return msa


def absolute2fractional(msa_model, slit, xposabs, yposabs):
    """
    Compute the fractional position in ``x`` within the slit in MSA coordinates.

    Parameters
    ----------
    msa_model : `~jwst.datamodels.MSAModel`
        NIRSpec MSA model.
    slit : `~jwst.datamodels.SlitModel`
        The slit data object.
    xposabs, yposabs : float
        (x, y) positions in the ``msa_frame``.

    Returns
    -------
    xpos : float
        The fractional X coordinates within the slit.
    """
    num, xcenter, ycenter, xsize, ysize = msa_model.Q5.data[slit.shutter_id]
    return (xposabs - xcenter) / (xsize / 2.)
