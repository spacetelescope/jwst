#
#  Module for 2d extraction
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import copy
import numpy as np
from astropy.io import fits
from astropy.modeling.models import Shift
from gwcs.utils import _toindex
from gwcs import wcstools

from .. import datamodels
from ..transforms import models as trmodels
from asdf import AsdfFile
from ..assign_wcs import nirspec


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract2d(input_model, which_subarray=None, apply_wavecorr=False, reffile=""):
    supported_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ', 'NRS_LAMP']
    
    exp_type = input_model.meta.exposure.type.upper()
    log.info('EXP_TYPE is {0}'.format(exp_type))

    wavecorr_supported_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ']
    if exp_type in wavecorr_supported_modes:
        if reffile and reffile.strip().upper() == 'N/A':
            apply_wavecorr = False
            log.info("WAVECORR reference file missing - skipping correction")
    else:
        apply_wavecorr = False
        log.info("Skipping wavecorr correction for EXP_TYPE {0}".format(exp_type))

    if exp_type not in supported_modes:
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model
    else:
        output_model = datamodels.MultiSlitModel()
        output_model.update(input_model)

    if exp_type in supported_modes:
        slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')
        # This is a cludge but will work for now.
        # This model keeps open_slits as an attribute.
        open_slits = slit2msa[1].slits[:]
        if which_subarray is not None:
            open_slits = [sub for sub in open_slits if sub.name==which_subarray]
        log.debug('open slits {0}'.format(open_slits))
        if apply_wavecorr and exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
            msa_model = get_msa_model(input_model)
        for slit in open_slits:
            slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
            xlo, xhi = _toindex(slit_wcs.bounding_box[0])
            ylo, yhi = _toindex(slit_wcs.bounding_box[1])

            # Add the slit offset to each slit WCS object
            tr = slit_wcs.get_transform('detector', 'sca')
            tr = Shift(xlo) & Shift(ylo) | tr
            slit_wcs.set_transform('detector', 'sca', tr.rename('dms2sca'))

            log.info('Name of subarray extracted: %s', slit.name)
            log.info('Subarray x-extents are: %s %s', xlo, xhi)
            log.info('Subarray y-extents are: %s %s', ylo, yhi)

            ext_data = input_model.data[ylo: yhi + 1, xlo: xhi + 1].copy()
            ext_err = input_model.err[ylo: yhi + 1, xlo: xhi + 1].copy()
            ext_dq = input_model.dq[ylo: yhi + 1, xlo: xhi + 1].copy()

            shape = ext_data.shape
            bounding_box= ((0, shape[1] - 1), (0, shape[0] - 1))
            slit_wcs.bounding_box = bounding_box

            # compute wavelengths
            x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
            ra, dec, lam = slit_wcs(x, y)
            if apply_wavecorr and _is_point_source(slit, exp_type, input_model.meta.target.source_type):
                if exp_type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
                    # pass lam = 2microns
                    source_xpos = get_source_xpos(input_model, slit, slit_wcs, lam=2, msa_model=msa_model)
                    aperture_name = slit.name
                else:
                    source_xpos = slit.source_xpos
                    # For the MSA the aperture name is "MOS"
                    aperture_name = "MOS"
                corr, dq_lam = compute_zero_point_correction(lam, reffile, source_xpos, aperture_name, slit_wcs)
                ## TODO: set a DQ flag to a TBD value for pixels where dq_lam == 0.
                ## The only purpose of dq_lam is to set that flag.
                lam = lam - corr
                log.info("Slit {0}: Wavelength zero-point correction applied.".format(slit.name))
            else:
                log.info("Slit {0}: Wavelength zero-point correction was not applied.".format(slit.name))
            new_model = datamodels.ImageModel(data=ext_data, err=ext_err, dq=ext_dq, wavelength=lam)
            new_model.meta.wcs = slit_wcs

            output_model.slits.append(new_model)
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            nslit = len(output_model.slits) - 1
            xlo_ind, xhi_ind, ylo_ind, yhi_ind = _toindex((xlo, xhi, ylo, yhi)).astype(np.int16)
            output_model.slits[nslit].name = str(slit.name)
            output_model.slits[nslit].xstart = xlo_ind + 1
            output_model.slits[nslit].xsize = (xhi_ind - xlo_ind) + 1
            output_model.slits[nslit].ystart = ylo_ind + 1
            output_model.slits[nslit].ysize = (yhi_ind - ylo_ind) + 1
            if exp_type.lower() == 'nrs_msaspec':
                output_model.slits[nslit].source_id = int(slit.source_id)
                output_model.slits[nslit].source_name = slit.source_name
                output_model.slits[nslit].source_alias = slit.source_alias
                output_model.slits[nslit].stellarity = float(slit.stellarity)
                output_model.slits[nslit].source_xpos = float(slit.source_xpos)
                output_model.slits[nslit].source_ypos = float(slit.source_ypos)
                output_model.slits[nslit].slitlet_id = int(slit.name)
                # for pathloss correction
                output_model.slits[nslit].nshutters = int(slit.nshutters)
    del input_model
    # Set the step status to COMPLETE
    output_model.meta.cal_step.extract_2d = 'COMPLETE'
    return output_model


# move to extract 2d for grism exposures
# add the wavelength range using the far left and far right orders?
# fselect1 = wrange_selector[0].index(input_model.meta.instrument.filter)
# fselect2 = wrange_selector[-1].index(input_model.meta.instrument.filter)
# lower_lam = wrange[0][fselect1][0]
# upper_lam = wrange[-1][fselect2][1]
# input_model.meta.wcsinfo.waverange_start = lower_lam
# input_model.meta.wcsinfo.waverange_end = upper_lam



def compute_zero_point_correction(lam, freference, source_xpos, aperture_name, slit_wcs):
    """
    Compute the Nirspec wavelength zero-point correction.
    
    Parameters
    ----------
    lam : nd-array like
        Wavelength array
    freference : str
        WAVECORR reference file name.
    source_xpos : float
        X position of the source as a fraction of the slit size.
    aperture_name : str
        Aperture name.
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

    Returns
    -------
    lambda_corr : ndarray like
        Wavelength correction.
    lam : ndarray like
        Interpolated wavelengths. Extrapolated values are reset to 0.
        This is returned so that the DQ array can be updated with a flag
        which indicates that no zero-point correction was done.
    """
    with datamodels.WaveCorrModel(freference) as wavecorr:
        for ap in wavecorr.apertures:
            if ap.aperture_name == aperture_name:
                log.info("Using wavelength zero-point correction for aperture {0}".format(ap.aperture_name))
                offset_model = ap.zero_point_offset.copy()
                variance = ap.variance.copy()
                width = ap.width
                break
        else:
            log.info("No wavelength zero-point correction found for slit {0}".format(aperture_name))
        
    dispersion = compute_dispersion(slit_wcs)
    deltax = source_xpos * width
    lam = lam.copy()
    l = lam[~np.isnan(lam)]
    offset_model.bounds_error = False
    correction = offset_model(l * 10 ** -6, [deltax]*l.size)
    lam[~np.isnan(lam)] = correction
    # The correction for pixels outside the slit and wavelengths
    # outside the wave_range is 0.
    lam[np.isnan(lam)] = 0.
    lambda_cor = dispersion * lam
    return lambda_cor, lam


def compute_dispersion(wcs):
    """
    Compute the pixel dispersion.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

    Returns
    -------
    dispersion : ndarray-like
        The pixel dispersion [in m].

    """
    xpix, ypix = wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1,1))
    xleft = xpix - 0.5
    xright = xpix + 0.5
    _, _, lamright = wcs(xright, ypix)
    _, _, lamleft = wcs(xleft, ypix)
    return (lamright - lamleft) * 10 ** -6


def _is_point_source(slit, exp_type, user_type):
    result = False
    if exp_type == 'NRS_MSASPEC':
        if slit.stellarity > 0.75:
            result = True
            log.info("Detected a point source in slit {0}, stelarity is {1}".format(slit.name, slit.stellarity))
        else:
            result = False
            log.info("Detected an extended source in slit {0}, stelarity is {1}".format(slit.name, slit.stellarity))
    else:
        # Get the value the user specified (if any)
        if (user_type is not None) and (user_type in ['POINT', 'EXTENDED']):
            # Use the value supplied by the user
            log.info('Detected a {0} source type in slit {1}'.format(user_type, slit.name))
            if user_type.strip().upper() == 'POINT':
                result = True
            else:
                result = False
        else:
            log.info("Unknown source type")

    return result

    
def get_source_xpos(input_model, slit, slit_wcs, lam, msa_model):
    """
    Compute the source position within the slit for a NIRSPEC FS.

    Parameters
    ----------
    input_model : `~jwst/datamodels/model_base.DataModel`
        The input to ``extract_2d``.
    slit : `~jwst/transforms/models/Slit`
        The slit tuple.
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    lam : float
        Wavelength in microns.

    Returns
    -------
    xpos : float
        X coordinate of the source as a fraction of the slit size.
    """
    xoffset = input_model.meta.dither.xoffset # in arcsec
    yoffset = input_model.meta.dither.yoffset # in arcsec
    v2ref = input_model.meta.wcsinfo.v2_ref # in arcsec
    v3ref = input_model.meta.wcsinfo.v3_ref # in arcsec
    v3idlyangle = input_model.meta.wcsinfo.v3yangle # in deg
    vparity = input_model.meta.wcsinfo.vparity
    
    idl2v23 = trmodels.IdealToV2V3(v3idlyangle, v2ref, v3ref, vparity)
    # Compute the location in V2,V3 [in arcsec]
    xv, yv = idl2v23(xoffset, yoffset)
    # The NIRSPEC transforms expect V2,V3 positions in deg
    xv = xv / 3600.
    yv = yv / 3600.
    
    v2v3_to_msa_frame = slit_wcs.get_transform("v2v3", "msa_frame")
    xpos_abs, ypos_abs, lam = v2v3_to_msa_frame(xv, yv, lam)
    xpos_frac = absolute2fractional(msa_model, slit, xpos_abs, ypos_abs)
    return xpos_frac


def get_msa_model(input_model):
    # Get the reference file used in constructing the WCS.
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
    input_model : `~jwst/datamodels/model_base.DataModel`
        The input to ``extract_2d``.
    slit : `~jwst/transforms/models/Slit`
        The slit tuple.
    xposabs, yposabs : float
        (x, y) positions in the ``msa_frame``.

    Returns
    -------
    xpos : float
        The fractional X coordinates within the slit.
    """
    num, xcenter, ycenter, xsize, ysize = msa_model.Q5.data[slit.shutter_id]
    return (xposabs - xcenter) / (xsize/ 2.)
