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

from jwst.assign_wcs import nirspec
from jwst.assign_wcs import util
from jwst.lib import pipe_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def nrs_extract2d(input_model, slit_names=None, source_ids=None):
    """
    Perform extract_2d calibration for NIRSpec exposures.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    slit_names : list containing strings or ints
        Slit names.
    source_ids : list containing strings or ints
        Source ids.

    Returns
    -------
    output_model : SlitModel or MultiSlitModel
        DataModel containing extracted slit(s)
    """
    exp_type = input_model.meta.exposure.type.upper()

    if (
        hasattr(input_model.meta.cal_step, "assign_wcs")
        and input_model.meta.cal_step.assign_wcs == "SKIPPED"
    ):
        log.info("assign_wcs was skipped")
        log.warning("extract_2d will be SKIPPED")
        input_model.meta.cal_step.extract_2d = "SKIPPED"
        return input_model

    if not (hasattr(input_model.meta, "wcs") and input_model.meta.wcs is not None):
        raise AttributeError(
            "Input model does not have a WCS object; assign_wcs should be run before extract_2d."
        )

    # Get all open slits from the WCS transforms
    open_slits = input_model.meta.wcs.get_transform("gwa", "slit_frame").slits

    # Select the slits to process
    open_slits = select_slits(open_slits, slit_names, source_ids)

    # NIRSpec BRIGHTOBJ (S1600A1 TSO) mode
    if exp_type == "NRS_BRIGHTOBJ":
        # the output model is a single SlitModel
        slit = open_slits[0]
        output_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit)
        set_slit_attributes(output_model, slit, xlo, xhi, ylo, yhi)
        try:
            get_source_xpos(output_model)
        except DitherMetadataError as e:
            log.warning(str(e))
            log.warning("Setting source position in slit to 0.0, 0.0")
            output_model.source_ypos = 0.0
            output_model.source_xpos = 0.0
        if "world" in input_model.meta.wcs.available_frames:
            orig_s_region = output_model.meta.wcsinfo.s_region.strip()
            util.update_s_region_nrs_slit(output_model)
            if orig_s_region != output_model.meta.wcsinfo.s_region.strip():
                log.info(f"extract_2d updated S_REGION to {output_model.meta.wcsinfo.s_region}")
    else:
        output_model = datamodels.MultiSlitModel()
        output_model.update(input_model)
        slits = []

        # Loop over all slit instances that are present
        for slit in open_slits:
            new_model, xlo, xhi, ylo, yhi = process_slit(input_model, slit)

            slits.append(new_model)
            orig_s_region = new_model.meta.wcsinfo.s_region.strip()
            # set x/ystart values relative to the image (screen) frame.
            # The overall subarray offset is recorded in model.meta.subarray.
            set_slit_attributes(new_model, slit, xlo, xhi, ylo, yhi)

            if new_model.meta.exposure.type.lower() == "nrs_fixedslit":
                if slit.name == input_model.meta.instrument.fixed_slit:
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
            if "world" in input_model.meta.wcs.available_frames:
                util.update_s_region_nrs_slit(new_model)
                if orig_s_region != new_model.meta.wcsinfo.s_region.strip():
                    log.info(f"Updated S_REGION to {new_model.meta.wcsinfo.s_region}")

            # Copy BUNIT values to output slit
            new_model.meta.bunit_data = input_model.meta.bunit_data
            new_model.meta.bunit_err = input_model.meta.bunit_err

        output_model.slits.extend(slits)

    return output_model


def select_slits(open_slits, slit_names, source_ids):
    """
    Select the slits to process.

    Parameters
    ----------
    open_slits : list
        List of open slits
    slit_names : list
        List of slit names to process
    source_ids : list
        List of source ids to process

    Returns
    -------
    selected_open_slits : list
        List of slits selected by slit_name or source_id
    """
    open_slit_names = [str(x.name) for x in open_slits]
    open_slit_source_ids = [str(x.source_id) for x in open_slits]
    selected_open_slits = []
    if slit_names is not None:
        matched_slits = []
        for this_slit in [str(x) for x in slit_names]:
            if this_slit in open_slit_names:
                matched_slits.append(this_slit)
            else:
                log.warn(f"Slit {this_slit} is not in the list of open slits.")
        for sub in open_slits:
            if str(sub.name) in matched_slits:
                selected_open_slits.append(sub)
    if source_ids is not None:
        matched_sources = []
        for this_id in [str(x) for x in source_ids]:
            if this_id in open_slit_source_ids:
                matched_sources.append(this_id)
            else:
                log.warn(f"Source id {this_id} is not in the list of open slits.")
        for sub in open_slits:
            if str(sub.source_id) in matched_sources:
                if sub not in selected_open_slits:
                    selected_open_slits.append(sub)
                else:
                    log.info(f"Source_id {sub.source_id} already selected (name {sub.name})")
    if len(selected_open_slits) > 0:
        log.info("Slits selected:")
        for this_slit in selected_open_slits:
            log.info(f"Name: {this_slit.name}, source_id: {this_slit.source_id}")
        return selected_open_slits
    else:
        log.info("All slits selected")
        return open_slits


def process_slit(input_model, slit):
    """
    Construct a data model for each slit.

    Extract the data.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model. The ``CubeModel`` is used only for TSO data, i.e.
        ``NRS_BRIGHTOBJ`` exposure or internal lamp exposures with lamp_mode
        set to ``BRIGHTOBJ``.
    slit : `~stdatamodels.jwst.transforms.models.Slit`
        A slit object.

    Returns
    -------
    new_model : `~jwst.datamodels.SlitModel`
        The new data model for a slit.
    xlo, xhi, ylo, yhi : float
        The corners of the extracted slit in pixel space.
    """
    new_model, xlo, xhi, ylo, yhi = extract_slit(input_model, slit)

    # Copy the DISPAXIS keyword to the output slit.
    new_model.meta.wcsinfo.dispersion_direction = input_model.meta.wcsinfo.dispersion_direction

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
    output_model.xsize = xhi - xlo
    output_model.ystart = ylo + 1  # account for FITS 1-indexed origin
    output_model.ysize = yhi - ylo
    output_model.source_id = int(slit.source_id)
    output_model.slit_ymin = slit.ymin
    output_model.slit_ymax = slit.ymax
    output_model.shutter_id = int(slit.shutter_id)  # for use in wavecorr
    log.debug(f"slit.ymin {slit.ymin}")
    if (
        output_model.meta.exposure.type.lower() in ["nrs_msaspec", "nrs_autoflat"]
        or output_model.meta.instrument.lamp_mode.upper == "MSASPEC"
    ):
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
        output_model.slit_xscale = float(slit.slit_xscale)
        output_model.slit_yscale = float(slit.slit_yscale)
        # for pathloss correction
        output_model.shutter_state = slit.shutter_state
    log.info("set slit_attributes completed")


def offset_wcs(slit_wcs):
    """
    Prepend a Shift transform to the slit WCS to account for subarrays.

    Parameters
    ----------
    slit_wcs : `~gwcs.wcs.WCS`
        The WCS for this  slit.

    Returns
    -------
    xlo, xhi, ylo, yhi : tuple of floats
        Indices of the bounding box of the WCS
    """
    xlo, xhi = _toindex(slit_wcs.bounding_box[0])
    ylo, yhi = _toindex(slit_wcs.bounding_box[1])

    # Add the slit offset to each slit WCS object
    tr = slit_wcs.get_transform("detector", "sca")
    tr = Shift(xlo) & Shift(ylo) | tr
    slit_wcs.set_transform("detector", "sca", tr.rename("dms2sca"))

    return xlo, xhi, ylo, yhi


def extract_slit(input_model, slit):
    """
    Extract a slit from a full frame image.

    Parameters
    ----------
    input_model : `~jwst.datamodels.image.ImageModel` or `~jwst.datamodels.cube.CubeModel`
        The input model.
    slit : `~stdatamodels.jwst.transforms.models.Slit`
        A slit object.

    Returns
    -------
    new_model : `~jwst.datamodels.SlitModel`
        The slit data model with WCS attached to it.
    """
    slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
    xlo, xhi, ylo, yhi = offset_wcs(slit_wcs)
    log.info(f"Name of subarray extracted: {slit.name}")
    log.info(f"Subarray x-extents are: {xlo} {xhi}")
    log.info(f"Subarray y-extents are: {ylo} {yhi}")
    ndim = len(input_model.data.shape)
    if ndim == 2:
        slit_slice = np.s_[ylo:yhi, xlo:xhi]
        ext_data = input_model.data[slit_slice].copy()
        ext_err = input_model.err[slit_slice].copy()
        ext_dq = input_model.dq[slit_slice].copy()
        ext_var_rnoise = input_model.var_rnoise[slit_slice].copy()
        ext_var_poisson = input_model.var_poisson[slit_slice].copy()
        int_times = None
    elif ndim == 3:
        slit_slice = np.s_[:, ylo:yhi, xlo:xhi]
        ext_data = input_model.data[slit_slice].copy()
        ext_err = input_model.err[slit_slice].copy()
        ext_dq = input_model.dq[slit_slice].copy()
        ext_var_rnoise = input_model.var_rnoise[slit_slice].copy()
        ext_var_poisson = input_model.var_poisson[slit_slice].copy()
        if pipe_utils.is_tso(input_model) and hasattr(input_model, "int_times"):
            log.debug("TSO data, so copying the INT_TIMES table.")
            int_times = input_model.int_times.copy()
        else:
            int_times = None
    else:
        raise ValueError(f"extract_2d does not work with {ndim} dimensional data")

    slit_wcs.bounding_box = util.wcs_bbox_from_shape(ext_data.shape)

    # compute wavelengths
    x, y = wcstools.grid_from_bounding_box(slit_wcs.bounding_box, step=(1, 1))
    ra, dec, lam = slit_wcs(x, y)
    lam = lam.astype(np.float32)
    new_model = datamodels.SlitModel(
        data=ext_data,
        err=ext_err,
        dq=ext_dq,
        wavelength=lam,
        var_rnoise=ext_var_rnoise,
        var_poisson=ext_var_poisson,
        int_times=int_times,
    )
    log.debug(f"Input model type is {str(input_model)}")
    new_model.update(input_model)
    new_model.meta.wcs = slit_wcs

    return new_model, xlo, xhi, ylo, yhi


class DitherMetadataError(Exception):
    """
    Handle DitherMetadataError exception.

    This exception is raised if a Slit object doesn't have the required
    dither attribute, or the x_ and y_ offsets aren't numeric.
    """

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
        raise DitherMetadataError(
            "meta.dither is not populated for the primary slit; "
            "Failed to estimate source position in slit."
        )

    if slit.meta.dither.x_offset is None or slit.meta.dither.y_offset is None:
        raise DitherMetadataError(
            "meta.dither.x(y)_offset values are None for primary slit; "
            "Failed to estimate source position in slit."
        )

    xoffset = slit.meta.dither.x_offset  # in arcsec
    yoffset = slit.meta.dither.y_offset  # in arcsec
    v2ref = slit.meta.wcsinfo.v2_ref  # in arcsec
    v3ref = slit.meta.wcsinfo.v3_ref  # in arcsec
    v3idlyangle = slit.meta.wcsinfo.v3yangle  # in deg
    vparity = slit.meta.wcsinfo.vparity

    idl2v23 = trmodels.IdealToV2V3(v3idlyangle, v2ref, v3ref, vparity)
    log.debug(f"wcsinfo: {v2ref}, {v3ref}, {v3idlyangle}, {vparity}")
    # Compute the location in V2,V3 [in arcsec]
    xv, yv = idl2v23(xoffset, yoffset)
    log.info(f"xoffset, yoffset, {xoffset}, {yoffset}")

    # Position in the virtual slit
    wavelength = 2.0  # microns, but it doesn't make any difference here
    xpos_slit, ypos_slit, lam_slit = slit.meta.wcs.get_transform("v2v3", "slit_frame")(
        xv, yv, wavelength
    )
    # Update slit.source_xpos, slit.source_ypos
    slit.source_xpos = xpos_slit
    slit.source_ypos = ypos_slit
    log.debug(f"Source X/Y position in V2V3: {xv}, {yv}")
    log.info(f"Source X/Y position in the slit: {xpos_slit}, {ypos_slit}")

    return xpos_slit
