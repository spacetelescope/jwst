import importlib
import logging

from gwcs.wcs import WCS

from jwst.assign_wcs.miri import store_dithered_position
from jwst.assign_wcs.util import (
    update_s_region_imaging,
    update_s_region_lrs,
    update_s_region_mrs,
    update_s_region_nrs_ifu,
    update_s_region_spectral,
)
from jwst.lib.dispaxis import get_dispersion_direction
from jwst.lib.exposure_types import IMAGING_TYPES, NRS_LAMP_MODE_SPEC_TYPES, SPEC_TYPES
from jwst.lib.wcs_utils import get_wavelengths

log = logging.getLogger(__name__)

__all__ = ["load_wcs"]


def load_wcs(input_model, reference_files=None, nrs_slit_y_range=None, nrs_ifu_slice_wcs=False):
    """
    Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The input data model. Updated in place.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
    nrs_slit_y_range : list
        The slit y-range for a NIRSpec slit. The center is (0, 0).
    nrs_ifu_slice_wcs : bool
        If True and the exposure type is NIRSpec IFU, then a full slice-based
        WCS that propagates slice IDs is produced.  This is intended primarily for
        diagnostic purposes.  If False and the exposure type is NIRSpec IFU,
        a slice map is internally applied to produce a fully coordinate-based
        WCS pipeline that does not require slice IDs on input.

    Returns
    -------
    input_model : `~stdatamodels.jwst.datamodels.JwstDataModel`
        The data model with the WCS information in the meta attribute.
    """
    if reference_files is not None:
        for ref_type, ref_file in reference_files.items():
            if ref_file not in ["N/A", ""]:
                reference_files[ref_type] = ref_file
            else:
                reference_files[ref_type] = None
    if (reference_files is None) or (not any(reference_files.values())):
        log.critical("assign_wcs needs reference files to compute the WCS, none were passed")
        raise ValueError("assign_wcs needs reference files to compute the WCS, none were passed")
    instrument = input_model.meta.instrument.name.lower()
    mod = importlib.import_module("." + instrument, "jwst.assign_wcs")

    if (
        input_model.meta.exposure.type.lower() in SPEC_TYPES
        or input_model.meta.instrument.lamp_mode.lower() in NRS_LAMP_MODE_SPEC_TYPES
    ):
        input_model.meta.wcsinfo.specsys = "BARYCENT"
        input_model.meta.wcsinfo.dispersion_direction = get_dispersion_direction(
            input_model.meta.exposure.type,
            input_model.meta.instrument.grating,
            input_model.meta.instrument.filter,
            input_model.meta.instrument.pupil,
        )
    if instrument.lower() == "nirspec":
        pipeline = mod.create_pipeline(input_model, reference_files, slit_y_range=nrs_slit_y_range)
    else:
        pipeline = mod.create_pipeline(input_model, reference_files)
    # Initialize the output model as a copy of the input
    # Make the copy after the WCS pipeline is created in order to pass updates to the model.
    if pipeline is None:
        input_model.meta.cal_step.assign_wcs = "SKIPPED"
        log.warning("assign_wcs: SKIPPED")
        return input_model

    wcs = WCS(pipeline)
    input_model.meta.wcs = wcs

    if (
        instrument.lower() == "nirspec"
        and input_model.meta.exposure.type.lower() not in IMAGING_TYPES
        and input_model.meta.instrument.grating.lower() != "mirror"
    ):
        cbbox = mod.generate_compound_bbox(input_model)
        input_model.meta.wcs.bounding_box = cbbox
    input_model.meta.cal_step.assign_wcs = "COMPLETE"
    exclude_types = [
        "nrc_wfss",
        "nrc_tsgrism",
        "nis_wfss",
        "nrs_fixedslit",
        "nrs_msaspec",
        "nrs_autowave",
        "nrs_autoflat",
        "nrs_lamp",
        "nrs_brightobj",
        "nis_soss",
        "mir_wfss",
    ]

    if input_model.meta.exposure.type.lower() not in exclude_types:
        imaging_types = IMAGING_TYPES.copy()
        imaging_types.update(["mir_lrs-slitless"])
        imaging_lrs_types = ["mir_lrs-fixedslit"]
        if input_model.meta.exposure.type.lower() in imaging_lrs_types:
            # uses slits corners in V2, V3 that are read in from the
            # lrs specwcs reference file
            update_s_region_lrs(input_model, reference_files)
            input_model.wavelength = get_wavelengths(input_model)
        elif input_model.meta.exposure.type.lower() in imaging_types:
            try:
                update_s_region_imaging(input_model)
            except Exception as exc:
                log.error(
                    f"Unable to update S_REGION for type {input_model.meta.exposure.type}: {exc}"
                )
            else:
                log.info(f"assign_wcs updated S_REGION to {input_model.meta.wcsinfo.s_region}")
            if input_model.meta.exposure.type.lower() == "mir_lrs-slitless":
                input_model.wavelength = get_wavelengths(input_model)
        elif input_model.meta.exposure.type.lower() == "nrs_ifu":
            update_s_region_nrs_ifu(input_model)

            # Attach a slice map in the regions attribute.
            # Optionally, use it to further revise the output WCS pipeline.
            mod.apply_slicemap(input_model, replace_wcs=(not nrs_ifu_slice_wcs))
        elif input_model.meta.exposure.type.lower() == "mir_mrs":
            update_s_region_mrs(input_model)
        else:
            try:
                update_s_region_spectral(input_model)
            except Exception as exc:
                log.info(
                    f"Unable to update S_REGION for type {input_model.meta.exposure.type}: {exc}"
                )

    # Store position of dithered pointing location in metadata for later spectral extraction
    if input_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        store_dithered_position(input_model)
        log.debug(
            "Storing dithered pointing location information: "
            f"{input_model.meta.dither.dithered_ra} {input_model.meta.dither.dithered_dec}"
        )
    log.info("COMPLETED assign_wcs")
    return input_model
