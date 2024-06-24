#
#  Module for 2d extraction of Nirspec fixed slits or MOS slitlets.
#
import logging
import numpy as np
from astropy.modeling.models import Shift
from gwcs.utils import _toindex
from gwcs import wcstools

from stdatamodels.jwst import datamodels
from stdatamodels.jwst.transforms import models as trmodels

from ..assign_wcs import nirspec
from ..assign_wcs import util
from ..lib import pipe_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def nrs_extract2d(input_model, slit_name=None):
    """
    Main extract_2d function for NIRSpec exposures.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    slit_name : str or int
        Slit name.
    """
    exp_type = input_model.meta.exposure.type.upper()

    if hasattr(input_model.meta.cal_step, 'assign_wcs') and input_model.meta.cal_step.assign_wcs == 'SKIPPED':
        log.info("assign_wcs was skipped")
        log.warning("extract_2d will be SKIPPED")
        input_model.meta.cal_step.extract_2d = "SKIPPED"
        return input_model

    if not (hasattr(input_model.meta, 'wcs') and input_model.meta.wcs is not None):
        raise AttributeError("Input model does not have a WCS object; assign_wcs should "
                             "be run before extract_2d.")

    slit2msa = input_model.meta.wcs.get_transform('slit_frame', 'msa_frame')
    # This is a kludge but will work for now.
    # This model keeps open_slits as an attribute.
    open_slits = slit2msa.slits[:]
    if slit_name is not None:
        new_open_slits = []
        slit_name = str(slit_name)
        for sub in open_slits:
            if str(sub.name) == slit_name:
                new_open_slits.append(sub)
                break
        if len(new_open_slits) == 0:
            raise AttributeError("Slit {} not in open slits.".format(slit_name))
        open_slits = new_open_slits

    # NIRSpec BRIGHTOBJ (S1600A1 TSO) mode
    if exp_type == 'NRS_BRIGHTOBJ':
        # the output model is a single SlitModel
        slit = open_slits[0]
        output_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit, exp_type)
        set_slit_attributes(output_model, slit, xlo, xhi, ylo, yhi)
        try:
            get_source_xpos(output_model)
        except DitherMetadataError as e:
            log.warning(str(e))
            log.warning("Setting source position in slit to 0.0, 0.0")
            output_model.source_ypos = 0.0
            output_model.source_xpos = 0.0
        if 'world' in input_model.meta.wcs.available_frames:
            orig_s_region = output_model.meta.wcsinfo.s_region.strip()
            util.update_s_region_nrs_slit(output_model)
            if orig_s_region != output_model.meta.wcsinfo.s_region.strip():
                log.info('extract_2d updated S_REGION to '
                         '{0}'.format(output_model.meta.wcsinfo.s_region))
    else:
        output_model = datamodels.MultiSlitModel()
        output_model.update(input_model)
        slits = []

        # Loop over all slit instances that are present
        for slit in open_slits:
            new_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit, exp_type)

            slits.append(new_model)
            orig_s_region = new_model.meta.wcsinfo.s_region.strip()
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            set_slit_attributes(new_model, slit, xlo, xhi, ylo, yhi)

            if (new_model.meta.exposure.type.lower() == 'nrs_fixedslit'):
                if (slit.name == input_model.meta.instrument.fixed_slit):
                    try:
                        get_source_xpos(new_model)
                    except DitherMetadataError as e:
                        log.warning(str(e))
                        log.warning("Setting source position in slit to 0.0, 0.0")
                        new_model.source_ypos = 0.0
                        new_model.source_xpos = 0.0
                else:
                    # ensure nonsense data never end up in non-primary slits
                    new_model.source_ypos = 0.0
                    new_model.source_xpos = 0.0

            # Update the S_REGION keyword value for the extracted slit
            if 'world' in input_model.meta.wcs.available_frames:
                util.update_s_region_nrs_slit(new_model)
                if orig_s_region != new_model.meta.wcsinfo.s_region.strip():
                    log.info(f'Updated S_REGION to {new_model.meta.wcsinfo.s_region}')

            # Copy BUNIT values to output slit
            new_model.meta.bunit_data = input_model.meta.bunit_data
            new_model.meta.bunit_err = input_model.meta.bunit_err

        output_model.slits.extend(slits)

    return output_model


def process_slit(input_model, slit, exp_type):
    """
    Construct a data model for each slit.

    Extract the data.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model. The ``CubeModel`` is used only for TSO data, i.e.
        ``NRS_BRIGHTOBJ`` exposure or internal lamp exposures with lamp_mode
        set to ``BRIGHTOBJ``.
    slit_name : str or int
        Slit name.
    exp_type : str
        The type of exposure. Supported types are
        ``NRS_FIXEDSLIT``, ``NRS_MSASPEC``, ``NRS_BRIGHTOBJ``, ``NRS_LAMP``,
        ``NRS_AUTOFLAT``, ``NRS_AUTOWAVE``

    Returns
    -------
    new_model : `~jwst.datamodels.SlitModel`
        The new data model for a slit.
    xlo, xhi, ylo, yhi : float
        The corners of the extracted slit in pixel space.

    """
    new_model, xlo, xhi, ylo, yhi = extract_slit(input_model, slit, exp_type)

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
        A `~stdatamodels.jwst.transforms.models.Slit` object representing a slit.
    xlo, xhi, ylo, yhi : float
        Indices into the data array where extraction should be done.
        These are converted to "pixel indices" - the center of a pixel.
    """
    output_model.name = str(slit.name)
    output_model.xstart = xlo + 1  # account for FITS 1-indexed origin
    output_model.xsize = (xhi - xlo)
    output_model.ystart = ylo + 1  # account for FITS 1-indexed origin
    output_model.ysize = (yhi - ylo)
    output_model.source_id = int(slit.source_id)
    output_model.slit_ymin = slit.ymin
    output_model.slit_ymax = slit.ymax
    output_model.shutter_id = int(slit.shutter_id)  # for use in wavecorr
    log.debug('slit.ymin {}'.format(slit.ymin))
    if output_model.meta.exposure.type.lower() in ['nrs_msaspec', 'nrs_autoflat'] or \
       output_model.meta.instrument.lamp_mode.upper == 'MSASPEC':
        # output_model.source_id = int(slit.source_id)
        output_model.source_name = slit.source_name
        output_model.source_alias = slit.source_alias
        output_model.stellarity = float(slit.stellarity)
        output_model.source_xpos = float(slit.source_xpos)
        output_model.source_ypos = float(slit.source_ypos)
        try:
            output_model.slitlet_id = int(slit.name)
        except ValueError:
            # Fixed slits in MOS data have string values for the name;
            # use the shutter ID instead
            output_model.slitlet_id = slit.shutter_id
        output_model.quadrant = int(slit.quadrant)
        output_model.xcen = int(slit.xcen)
        output_model.ycen = int(slit.ycen)
        output_model.dither_position = int(slit.dither_position)
        output_model.source_ra = float(slit.source_ra)
        output_model.source_dec = float(slit.source_dec)
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
    slit : `~stdatamodels.jwst.transforms.models.Slit`
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
    log.info(f'Name of subarray extracted: {slit.name}')
    log.info(f'Subarray x-extents are: {xlo} {xhi}')
    log.info(f'Subarray y-extents are: {ylo} {yhi}')
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
        if pipe_utils.is_tso(input_model) and hasattr(input_model, 'int_times'):
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
    log.debug(f'Input model type is {input_model.__class__.__name__}')
    new_model.update(input_model)
    new_model.meta.wcs = slit_wcs

    return new_model, xlo, xhi, ylo, yhi


class DitherMetadataError(Exception):
    pass


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

    if not hasattr(slit.meta, "dither"):
        raise DitherMetadataError('meta.dither is not populated for the primary slit; '
                                  'Failed to estimate source position in slit.')

    if slit.meta.dither.x_offset is None or slit.meta.dither.y_offset is None:
        raise DitherMetadataError('meta.dither.x(y)_offset values are None for primary slit; '
                                  'Failed to estimate source position in slit.')

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
    wavelength = 2.0 # microns, but it doesn't make any difference here
    xpos_slit, ypos_slit, lam_slit = slit.meta.wcs.get_transform('v2v3', 'slit_frame')(
        xv, yv, wavelength)
    # Update slit.source_xpos, slit.source_ypos
    slit.source_xpos = xpos_slit
    slit.source_ypos = ypos_slit
    log.debug('Source X/Y position in V2V3: {0}, {1}'.format(xv, yv))
    log.info('Source X/Y position in the slit: {0}, {1}'.format(xpos_slit, ypos_slit))

    return xpos_slit
