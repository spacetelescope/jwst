"""
Module for applying wavelength corrections to NIRSpec MOS and FS
slits in which the source is off center. The correction is applied
only to point sources.

The algorithm uses a reference file which is a look up table of
wavelength_correction as a function of slit_x_position and wavelength.
The ``x`` direction is the one parallel to dispersion/wavelength for
both MOS and FS slits.

The slit_x_position is determined in the following way.

MOS data:
``slit.source_xpos`` is used. It is read in from the msa_metadata_file
in the assign_wcs step.

FS data:
``input_model.meta.dither.x_offset``, populated by SDP, is used. The
value is in the telescope Ideal frame. To calculate the slit position,
it is transformed Ideal --> V2, V3 --> slit_frame.

The interpolation of the lookup table uses the absolute fractional offset in x.
This can be computed in two ways. The above transform gives the absolute
fractional position within the slit directly. The second way is to transform the
Ideal source_xpos position to the MSA frame and use the metrology data in the
``msa`` reference file to compute ``(xposabs - xcenter) / (xsize)``.
Both ways give the same result, the current code uses the first one.

"""
import logging
import numpy as np
from gwcs import wcstools
from gwcs import coordinate_frames as cf
from astropy import units as u
from astropy.modeling import tabular
from astropy.modeling.mappings import Identity

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms import models as trmodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def do_correction(input_model, wavecorr_file):
    """ Main wavecorr correction function for NIRSpec MOS and FS.

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

    # For BRIGHTOBJ, operate on the single SlitModel
    if isinstance(input_model, datamodels.SlitModel):
        if _is_point_source(input_model, exp_type):
            apply_zero_point_correction(output_model, wavecorr_file)
    else:
        # For FS only work on point source slits with
        # position information
        corrected = False
        if exp_type == 'NRS_FIXEDSLIT':
            primary_slit = input_model.meta.instrument.fixed_slit
            for slit in output_model.slits:
                if _is_point_source(slit, exp_type):
                    # If fixed slit was not defined via MSA file,
                    # it must have dither information to find the
                    # source position, and it must be the primary slit
                    if not _is_msa_fixed_slit(slit):
                        if slit.name != primary_slit:
                            log.warning(f'Skipping wavecorr correction for '
                                        f'non-primary slit {slit.name}')
                            continue
                        if not hasattr(slit.meta, "dither"):
                            log.warning('meta.dither is not populated for the primary slit')
                            log.warning('Skipping wavecorr correction')
                            continue
                        if slit.meta.dither.x_offset is None or slit.meta.dither.y_offset is None:
                            log.warning('dither.x(y)_offset values are None for primary slit')
                            log.warning('Skipping wavecorr correction')
                            input_model.meta.cal_step.wavecorr = 'SKIPPED'
                            continue
                    completed = apply_zero_point_correction(slit, wavecorr_file)
                    if completed:
                        corrected = True
                        slit.meta.cal_step.wavecorr = 'COMPLETE'
                    else:  # pragma: no cover
                        log.warning(f'Corrections are not invertible for slit {slit.name}')
                        log.warning('Skipping wavecorr correction')
                        slit.meta.cal_step.wavecorr = 'SKIPPED'

            if corrected:
                output_model.meta.cal_step.wavecorr = 'COMPLETE'
            else:
                output_model.meta.cal_step.wavecorr = 'SKIPPED'

        # For MOS work on all slits containing a point source
        else:
            for slit in output_model.slits:
                if _is_point_source(slit, exp_type):
                    completed = apply_zero_point_correction(slit, wavecorr_file)
                    if completed:
                        slit.meta.cal_step.wavecorr = 'COMPLETE'
                        corrected = True
                    else: # pragma: no cover
                        log.warning(f'Corrections are not invertible for slit {slit.name}')
                        log.warning('Skipping wavecorr correction')
                        slit.meta.cal_step.wavecorr = 'SKIPPED'
                else:
                    slit.meta.cal_step.wavecorr = 'SKIPPED'

            if corrected:
                output_model.meta.cal_step.wavecorr = 'COMPLETE'
            else:
                output_model.meta.cal_step.wavecorr = 'SKIPPED'

    return output_model


def apply_zero_point_correction(slit, reffile):
    """ Apply the NIRSpec wavelength zero-point correction.

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
    log.info(f'slit name {slit.name}')
    slit_wcs = slit.meta.wcs

    # Get the source position in the slit and set the aperture name
    if slit.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        # Check for fixed slits defined via MSA files in
        # MOS/FS combination processing: they should not have their
        # source position overridden by dither keywords
        if _is_msa_fixed_slit(slit):
            # get the planned source position
            source_xpos = slit.source_xpos
        else:
            # get the source position from the dither offsets
            source_xpos = get_source_xpos(slit)
        aperture_name = slit.name
    else:
        source_xpos = slit.source_xpos
        # For the MSA the aperture name is "MOS"
        aperture_name = "MOS"

    lam = slit.wavelength.copy() * 1e-6
    dispersion = compute_dispersion(slit.meta.wcs)

    wave2wavecorr = calculate_wavelength_correction_transform(
        lam, dispersion, reffile, source_xpos, aperture_name)

    # wave2wavecorr should not be None for real data
    if wave2wavecorr is None: # pragma: no cover
        completed = False
        return completed
    else:        
        # Make a new frame to insert into the slit wcs object
        slit_spatial = cf.Frame2D(name='slit_spatial', axes_order=(0, 1), 
                                  unit=("", ""), axes_names=('x_slit', 'y_slit'))
        spec = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=(u.micron,),
                                axes_names=('wavelength',))
        wcorr_frame = cf.CompositeFrame(
            [slit_spatial, spec], name='wavecorr_frame')
        
        # Insert the new transform into the slit wcs object
        wave2wavecorr = Identity(2) & wave2wavecorr
        slit_wcs.insert_frame('slit_frame', wave2wavecorr, wcorr_frame)
        
        # Update the stored wavelengths for the slit
        slit.wavelength = compute_wavelength(slit_wcs)
        
        completed = True
        return completed
    

def calculate_wavelength_correction_transform(lam, dispersion, freference, 
                                              source_xpos, aperture_name):
    """ Generate a WCS transform for the NIRSpec wavelength zero-point correction
    and add it to the WCS for each slit.

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
                log.info(f'Using wavelength zero-point correction for aperture {ap.aperture_name}')
                offset_model = ap.zero_point_offset.copy()
                break
        else:
            log.info(f'No wavelength zero-point correction found for slit {aperture_name}')
        
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
        wave2wavecorr = tabular.Tabular1D(points=lam_mean, 
                                          lookup_table=lam_corrected, 
                                          bounds_error=False, 
                                          fill_value=None, 
                                          name='wave2wavecorr')
    
        return wave2wavecorr
    
    else:
        # output wavelengths are not monotonically increasing
        return None


def compute_dispersion(wcs, xpix=None, ypix=None):
    """ Compute the pixel dispersion.

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


def compute_wavelength(wcs, xpix=None, ypix=None):
    """ Compute the pixel wavelength.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    xpix : ndarray, float, optional
    ypix : ndarray, float, optional
        Compute the wavelength at the x, y pixels.
        If not provided the dispersion is computed on a
        grid based on ``wcs.bounding_box``.

    Returns
    -------
    wavelength : ndarray
        The wavelength [in microns].

    """
    if xpix is None or ypix is None:
        xpix, ypix = wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1, 1))
        
    _, _, lam = wcs(xpix, ypix)
    return lam


def _is_msa_fixed_slit(slit):
    """
    Determine if a fixed slit source was defined via a MSA file.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.transforms.models.Slit`
        A slit object.
    """
    # Fixed slits defined via MSA files in  MOS/FS combination
    # processing will have a non-empty source name
    if (not hasattr(slit, 'source_name')
            or slit.source_name is None
            or slit.source_name == ""):
        return False
    else:
        return True


def _is_point_source(slit, exp_type):
    """
    Determine if a source is a point source.

    Parameters
    ----------
    slit : `~stdatamodels.jwst.transforms.models.Slit`
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

    if src_type is not None and src_type.upper() in ['POINT', 'EXTENDED']:
        # Use the supplied value
        log.info(f'Detected a {src_type} source type in slit {slit.name}')
        if src_type.strip().upper() == 'POINT':
            result = True
        else:
            result = False
    else:
        log.info("Unknown source type")

    return result


def get_source_xpos(slit):
    """
    Compute the source position within the slit for a NIRSpec fixed slit.

    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`
        The slit object.

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
    log.debug("wcsinfo: {0}, {1}, {2}, {3}".format(v2ref, v3ref, v3idlyangle, vparity))
    # Compute the location in V2,V3 [in arcsec]
    xv, yv = idl2v23(xoffset, yoffset)
    log.info(f'xoffset, yoffset, {xoffset}, {yoffset}')

    # Position in the virtual slit
    xpos_slit, ypos_slit, lam_slit = slit.meta.wcs.get_transform('v2v3', 'slit_frame')(
        xv, yv, 2)
    # Update slit.source_xpos, slit.source_ypos
    slit.source_xpos = xpos_slit
    slit.source_ypos = ypos_slit
    log.debug('Source X/Y position in V2V3: {0}, {1}'.format(xv, yv))
    log.info('Source X/Y position in the slit: {0}, {1}'.format(xpos_slit, ypos_slit))

    return xpos_slit
