#
#  Module for 2d extraction of Nirspec fixed slits or MOS slitlets.
#
import logging
import warnings
import numpy as np
from astropy.modeling.models import Shift
from gwcs.utils import _toindex
from gwcs import wcstools

from .. import datamodels
from ..transforms import models as trmodels
from ..assign_wcs import nirspec
from ..assign_wcs import util
from ..lib import pipe_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def nrs_extract2d(input_model, slit_name=None, apply_wavecorr=False, reference_files={}):
    """
    Main extract_2d function for Nirspec exposures.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    slit_name : str or int
        Slit name.
    apply_wavecorr : bool
        Flag whether to apply the zero point wavelength correction to
        Nirspec exposures.
    reference_files : dict
        Reference files - uses the ``wavecorr`` reference file.
    """
    exp_type = input_model.meta.exposure.type.upper()

    wavecorr_supported_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                                'NRS_AUTOFLAT']

    if exp_type in wavecorr_supported_modes:
        reffile = reference_files['wavecorr']
        if reffile and reffile.strip().upper() == 'N/A':
            apply_wavecorr = False
            warnings.warn("WAVECORR reference file missing - skipping correction")
    else:
        reffile = None
        apply_wavecorr = False
        log.info("Skipping wavecorr correction for EXP_TYPE {0}".format(exp_type))

    if hasattr(input_model.meta.cal_step, 'assign_wcs') and input_model.meta.cal_step.assign_wcs == 'SKIPPED':
        log.info("assign_wcs was skipped")
        log.warning("extract_2d: SKIPPED")
        input_model.meta.cal_step.extract_2d = "SKIPPED"
        return input_model

    if not (hasattr(input_model.meta, 'wcs') and input_model.meta.wcs is not None):
        raise AttributeError("Input model does not have a WCS object; assign_wcs should "
                             "be run before extract_2d.")

    slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')
    # This is a cludge but will work for now.
    # This model keeps open_slits as an attribute.
    open_slits = slit2msa.slits[:]
    if slit_name is not None:
        open_slits = [sub for sub in open_slits if sub.name == slit_name]
    log.debug('open slits {0}'.format(open_slits))
    if exp_type == 'NRS_BRIGHTOBJ':
        # the output model is a SlitModel
        slit = open_slits[0]
        output_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit,
                                                        exp_type, apply_wavecorr, reffile)
        set_slit_attributes(output_model, slit, xlo, xhi, ylo, yhi)
        orig_s_region = output_model.meta.wcsinfo.s_region.strip()
        util.update_s_region_spectral(output_model)
        if orig_s_region != output_model.meta.wcsinfo.s_region.strip():
            log.info('extract_2d updated S_REGION to '
                     '{0}'.format(output_model.meta.wcsinfo.s_region))
    else:
        output_model = datamodels.MultiSlitModel()
        output_model.update(input_model)
        slits = []
        for slit in open_slits:
            new_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit,
                                                         exp_type, apply_wavecorr, reffile)

            slits.append(new_model)
            orig_s_region = new_model.meta.wcsinfo.s_region.strip()
            util.update_s_region_spectral(new_model)
            if orig_s_region != new_model.meta.wcsinfo.s_region.strip():
                log.info('extract_2d updated S_REGION to {0}'.format(new_model.meta.wcsinfo.s_region))
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            set_slit_attributes(new_model, slit, xlo, xhi, ylo, yhi)

            # Copy BUNIT values to output slit
            new_model.meta.bunit_data = input_model.meta.bunit_data
            new_model.meta.bunit_err = input_model.meta.bunit_err
        output_model.slits.extend(slits)
    return output_model


def process_slit(input_model, slit, exp_type, apply_wavecorr, reffile):
    """
    Construct a data model for each slit.

    Extract the data and apply the wavelength
    zero point correction if requested.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model. The ``CubeModel`` is used only for TSO data, i.e.
        ``NRS_BRIGHTOBJ`` exposure.
    slit_name : str or int
        Slit name.
    exp_type : str
        The type of exposure. Supported types are
        ``NRS_FIXEDSLIT``, ``NRS_MSASPEC``, ``NRS_BRIGHTOBJ``
    apply_wavecorr : bool
        Flag whether to apply the zero point wavelength correction.
    reffile : str
        Path to ``wavecorr`` reference file.

    Returns
    -------
    new_model : `~jwst/datamodels/SlitModel`
        The new data model for a slit.
    xlo, xhi, ylo, yhi : float
        The corners of the extracted slit in pixel space.

    """
    new_model, xlo, xhi, ylo, yhi = extract_slit(input_model, slit, exp_type)
    if apply_wavecorr and _is_point_source(slit, exp_type, input_model.meta.target.source_type):
        apply_zero_point_correction(new_model, slit, reffile)
        log.info("Slit {0}: Wavelength zero-point correction applied.".format(slit.name))
    else:
        log.info("Slit {0}: Wavelength zero-point correction "
                 "was not applied.".format(slit.name))
    # Copy the DISPAXIS keyword to the output slit.
    new_model.meta.wcsinfo.dispersion_direction = \
        input_model.meta.wcsinfo.dispersion_direction
    return new_model, xlo, xhi, ylo, yhi


def set_slit_attributes(output_model, slit, xlo, xhi, ylo, yhi):
    """
    Set the slit attributes.

    Parameters
    ----------
    output_model : `~jwst.datamodels.SlitModel`
        The output model representing a slit.
    slit : namedtuple
        A `~jwst.transforms.models.Slit` object representing a slit.
    xlo, xhi, ylo, yhi : float
        Indices into the data array where extraction should be done.
        These are converted to "pixel indices" - the center of a pixel.
    """
    output_model.name = str(slit.name)
    output_model.xstart = xlo + 1 # FITS 1-indexed
    output_model.xsize = (xhi - xlo)
    output_model.ystart = ylo + 1 # FITS 1-indexed
    output_model.ysize = (yhi - ylo)
    output_model.source_id = int(slit.source_id)
    if output_model.meta.exposure.type.lower() in ['nrs_msaspec', 'nrs_autoflat']:
        #output_model.source_id = int(slit.source_id)
        output_model.source_name = slit.source_name
        output_model.source_alias = slit.source_alias
        output_model.stellarity = float(slit.stellarity)
        output_model.source_xpos = float(slit.source_xpos)
        output_model.source_ypos = float(slit.source_ypos)
        output_model.slitlet_id = int(slit.name)
        output_model.quadrant = int(slit.quadrant)
        output_model.xcen = int(slit.xcen)
        output_model.ycen = int(slit.ycen)
        output_model.dither_position = int(slit.dither_position)
        # for pathloss correction
        output_model.shutter_state = slit.shutter_state
    log.info('set slit_attributes completed')


def offset_wcs(slit_wcs):
    """
    Prepend a Shift transform to the slit WCS to account for subarrays.

    Parameters
    ----------
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS for this  slit.
    slit_name : str
        The name of the slit.
    """
    xlo, xhi = _toindex(slit_wcs.bounding_box[0])
    ylo, yhi = _toindex(slit_wcs.bounding_box[1])

    # Add the slit offset to each slit WCS object
    tr = slit_wcs.get_transform('detector', 'sca')
    tr = Shift(xlo) & Shift(ylo) | tr
    slit_wcs.set_transform('detector', 'sca', tr.rename('dms2sca'))

    return xlo, xhi, ylo, yhi


def extract_slit(input_model, slit, exp_type):
    """
    Extract a slit from a full frame image.

    Parameters
    ----------
    input_model : `~jwst.datamodels.image.ImageModel` or `~jwst.datamodels.cube.CubeModel`
        The input model.
    slit : `~jwst.transforms.models.Slit`
        A slit object.
    exp_type : str
        The exposure type.

    Returns
    -------
    new_model : `~jwst.datamodels.SlitModel`
        The slit data model with WCS attached to it.
    """
    slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
    xlo, xhi, ylo, yhi = offset_wcs(slit_wcs)
    log.info('Name of subarray extracted: %s', slit.name)
    log.info('Subarray x-extents are: %s %s', xlo, xhi)
    log.info('Subarray y-extents are: %s %s', ylo, yhi)
    ndim = len(input_model.data.shape)
    if ndim == 2:
        slit_slice = np.s_[ylo: yhi, xlo: xhi]
        ext_data = input_model.data[slit_slice].copy()
        ext_err = input_model.err[slit_slice].copy()
        ext_dq = input_model.dq[slit_slice].copy()
        ext_var_rnoise = input_model.var_rnoise[slit_slice].copy()
        ext_var_poisson = input_model.var_poisson[slit_slice].copy()
        int_times = None
    elif ndim == 3:
        slit_slice = np.s_[:, ylo: yhi, xlo: xhi]
        ext_data = input_model.data[slit_slice].copy()
        ext_err = input_model.err[slit_slice].copy()
        ext_dq = input_model.dq[slit_slice].copy()
        ext_var_rnoise = input_model.var_rnoise[slit_slice].copy()
        ext_var_poisson = input_model.var_poisson[slit_slice].copy()
        if (pipe_utils.is_tso(input_model) and
            hasattr(input_model, 'int_times')):
                log.debug("TSO data, so copying the INT_TIMES table.")
                int_times = input_model.int_times.copy()
        else:
                int_times = None
    else:
        raise ValueError("extract_2d does not work with "
                         "{0} dimensional data".format(ndim))

    slit_wcs.bounding_box = util.wcs_bbox_from_shape(ext_data.shape)

    # compute wavelengths
    x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
    ra, dec, lam = slit_wcs(x, y)
    lam = lam.astype(np.float32)
    new_model = datamodels.SlitModel(data=ext_data, err=ext_err, dq=ext_dq, wavelength=lam,
                                     var_rnoise=ext_var_rnoise, var_poisson=ext_var_poisson,
                                     int_times=int_times)
    log.info('Input model type is {}'.format(input_model.__class__.__name__))
    new_model.update(input_model)
    new_model.meta.wcs = slit_wcs
    return new_model, xlo, xhi, ylo, yhi


def apply_zero_point_correction(model, slit, reffile):
    """
    Apply the NIRSpec wavelength zero-point correction.

    Parameters
    ----------
    model : `~jwst.datamodels.SlitModel`, `~jwst.datamodels.cube.CubeModel`
        The output of `extract_slit`.
    slit : `~jwst.transforms.models.Slit`
        A slit object.
    reffile : str
        The ``wavecorr`` reference file.
    """
    slit_wcs = model.meta.wcs

    if model.meta.exposure.type in ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ']:
        # pass lam = 2 microns
        # needed for wavecorr with fixed slits
        msa_model = get_msa_model(model)
        source_xpos = get_source_xpos(model, slit, slit_wcs, lam=2,
                                      msa_model=msa_model)
        aperture_name = slit.name
    else:
        source_xpos = slit.source_xpos
        # For the MSA the aperture name is "MOS"
        aperture_name = "MOS"
    lam = model.wavelength
    dispersion = compute_dispersion(model.meta.wcs)
    corr, dq_lam = compute_zero_point_correction(lam, reffile, source_xpos,
                                                 aperture_name, dispersion)
    ## TODO: set a DQ flag to a TBD value for pixels where dq_lam == 0.
    ## The only purpose of dq_lam is to set that flag.
    model.wavelength = lam - corr


def compute_zero_point_correction(lam, freference, source_xpos, aperture_name, dispersion):
    """
    Compute the Nirspec wavelength zero-point correction.

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
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

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
                log.info("Using wavelength zero-point correction for aperture {0}".format(ap.aperture_name))
                offset_model = ap.zero_point_offset.copy()
                # TODO: implement variance
                #variance = ap.variance.copy()
                width = ap.width
                break
        else:
            log.info("No wavelength zero-point correction found for slit {0}".format(aperture_name))

    deltax = source_xpos * width
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


def compute_dispersion(wcs):
    """
    Compute the pixel dispersion.

    Parameters
    ----------
    wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.

    Returns
    -------
    dispersion : ndarray
        The pixel dispersion [in m].

    """
    xpix, ypix = wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1, 1))
    xleft = xpix - 0.5
    xright = xpix + 0.5
    _, _, lamright = wcs(xright, ypix)
    _, _, lamleft = wcs(xleft, ypix)
    return (lamright - lamleft) * 10 ** -6


def _is_point_source(slit, exp_type, user_type):
    """
    Determine if a source is a point source.

    Parameters
    ----------
    slit : `~jwst.transforms.models.Slit`
        A slit object.
    exp_type : str
        The exposure type
    user_type : str
        User defined source type (coming from the proposal).
    """
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
        if (user_type is not None) and (user_type.upper() in ['POINT', 'EXTENDED']):
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
    Compute the source position within the slit for a NIRSPEC fixed slit.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel`
        The input to ``extract_2d``.
    slit : `~jwst.transforms.models.Slit`
        The slit object (tuple).
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS object for this slit.
    lam : float
        Wavelength in microns.

    Returns
    -------
    xpos : float
        X coordinate of the source as a fraction of the slit size.
    """
    xoffset = input_model.meta.dither.x_offset # in arcsec
    yoffset = input_model.meta.dither.y_offset # in arcsec
    v2ref = input_model.meta.wcsinfo.v2_ref # in arcsec
    v3ref = input_model.meta.wcsinfo.v3_ref # in arcsec
    v3idlyangle = input_model.meta.wcsinfo.v3yangle # in deg
    vparity = input_model.meta.wcsinfo.vparity

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
    input_model : `~jwst.datamodels.ImageModel`
        The input to ``extract_2d``.
    slit : `~jwst.transforms.models.Slit`
        The slit object.
    xposabs, yposabs : float
        (x, y) positions in the ``msa_frame``.

    Returns
    -------
    xpos : float
        The fractional X coordinates within the slit.
    """
    num, xcenter, ycenter, xsize, ysize = msa_model.Q5.data[slit.shutter_id]
    return (xposabs - xcenter) / (xsize / 2.)
