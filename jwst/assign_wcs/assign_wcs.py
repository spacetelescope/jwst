import logging
import importlib
from gwcs.wcs import WCS
from .util import (
    update_s_region_spectral,
    update_s_region_imaging,
    update_s_region_nrs_ifu,
    update_s_region_mrs,
    update_s_region_lrs,
)
from jwst.lib.exposure_types import IMAGING_TYPES, SPEC_TYPES, NRS_LAMP_MODE_SPEC_TYPES
from jwst.lib.dispaxis import get_dispersion_direction
from jwst.lib.wcs_utils import get_wavelengths
from .miri import store_dithered_position

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["load_wcs"]


def load_wcs(input_model, reference_files=None, nrs_slit_y_range=None, nrs_ifu_slice_wcs=False):
    """
    Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        The input data model.
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
    output_model : `~jwst.datamodels.JwstDataModel`
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

    output_model = input_model.copy()
    wcs = WCS(pipeline)
    output_model.meta.wcs = wcs
    if (
        instrument.lower() == "nirspec"
        and output_model.meta.exposure.type.lower() not in IMAGING_TYPES
    ):
        cbbox = mod.generate_compound_bbox(output_model)
        output_model.meta.wcs.bounding_box = cbbox
    output_model.meta.cal_step.assign_wcs = "COMPLETE"
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
    ]

    if output_model.meta.exposure.type.lower() not in exclude_types:
        imaging_types = IMAGING_TYPES.copy()
        imaging_types.update(["mir_lrs-slitless"])
        imaging_lrs_types = ["mir_lrs-fixedslit"]
        if output_model.meta.exposure.type.lower() in imaging_lrs_types:
            # uses slits corners in V2, V3 that are read in from the
            # lrs specwcs reference file
            update_s_region_lrs(output_model, reference_files)
            output_model.wavelength = get_wavelengths(output_model)
        elif output_model.meta.exposure.type.lower() in imaging_types:
            try:
                update_s_region_imaging(output_model)
            except Exception as exc:
                log.error(
                    f"Unable to update S_REGION for type {output_model.meta.exposure.type}: {exc}"
                )
            else:
                log.info(f"assign_wcs updated S_REGION to {output_model.meta.wcsinfo.s_region}")
            if output_model.meta.exposure.type.lower() == "mir_lrs-slitless":
                output_model.wavelength = get_wavelengths(output_model)
        elif output_model.meta.exposure.type.lower() == "nrs_ifu":
            update_s_region_nrs_ifu(output_model)

            # Attach a slice map in the regions attribute.
            # Optionally, use it to further revise the output WCS pipeline.
            mod.apply_slicemap(output_model, replace_wcs=(not nrs_ifu_slice_wcs))
        elif output_model.meta.exposure.type.lower() == "mir_mrs":
            update_s_region_mrs(output_model)
        else:
            try:
                update_s_region_spectral(output_model)
            except Exception as exc:
                log.info(
                    f"Unable to update S_REGION for type {output_model.meta.exposure.type}: {exc}"
                )

    # Store position of dithered pointing location in metadata for later spectral extraction
    if output_model.meta.exposure.type.lower() == "mir_lrs-fixedslit":
        store_dithered_position(output_model)
        log.debug(
            "Storing dithered pointing location information: "
            f"{output_model.meta.dither.dithered_ra} {output_model.meta.dither.dithered_dec}"
        )
    log.info("COMPLETED assign_wcs")
    return output_model
