import logging
import importlib
from gwcs.wcs import WCS
from .util import (update_s_region_spectral, update_s_region_imaging,
                   update_s_region_nrs_ifu, update_s_region_mrs)
from ..lib.exposure_types import IMAGING_TYPES, SPEC_TYPES, NRS_LAMP_MODE_SPEC_TYPES
from ..lib.dispaxis import get_dispersion_direction
from ..lib.wcs_utils import get_wavelengths

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["load_wcs"]


def load_wcs(input_model, reference_files={}, nrs_slit_y_range=None):
    """
    Create a gWCS object and store it in ``Model.meta``.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        The exposure.
    reference_files : dict
        A dict {reftype: reference_file_name} containing all
        reference files that apply to this exposure.
    nrs_slit_y_range : list
        The slit y-range for a Nirspec slit. The center is (0, 0).
    """
    if reference_files:
        for ref_type, ref_file in reference_files.items():
            if ref_file not in ["N/A", ""]:
                reference_files[ref_type] = ref_file
            else:
                reference_files[ref_type] = None
    if not any(reference_files.values()):
        log.critical("assign_wcs needs reference files to compute the WCS, none were passed")
        raise ValueError("assign_wcs needs reference files to compute the WCS, none were passed")
    instrument = input_model.meta.instrument.name.lower()
    mod = importlib.import_module('.' + instrument, 'jwst.assign_wcs')

    if input_model.meta.exposure.type.lower() in SPEC_TYPES or \
       input_model.meta.instrument.lamp_mode.lower() in NRS_LAMP_MODE_SPEC_TYPES:
        input_model.meta.wcsinfo.specsys = "BARYCENT"
        input_model.meta.wcsinfo.dispersion_direction = \
            get_dispersion_direction(
                input_model.meta.exposure.type,
                input_model.meta.instrument.grating,
                input_model.meta.instrument.filter,
                input_model.meta.instrument.pupil)

    if instrument.lower() == 'nirspec':
        pipeline = mod.create_pipeline(input_model, reference_files, slit_y_range=nrs_slit_y_range)
    else:
        pipeline = mod.create_pipeline(input_model, reference_files)
    # Initialize the output model as a copy of the input
    # Make the copy after the WCS pipeline is created in order to pass updates to the model.
    if pipeline is None:
        input_model.meta.cal_step.assign_wcs = 'SKIPPED'
        log.warning("assign_wcs: SKIPPED")
        return input_model
    else:
        output_model = input_model.copy()
        wcs = WCS(pipeline)
        output_model.meta.wcs = wcs
        output_model.meta.cal_step.assign_wcs = 'COMPLETE'
        exclude_types = ['nrc_wfss', 'nrc_tsgrism', 'nis_wfss',
                         'nrs_fixedslit', 'nrs_msaspec',
                         'nrs_autowave', 'nrs_autoflat', 'nrs_lamp',
                         'nrs_brightobj', 'nis_soss']

        if output_model.meta.exposure.type.lower() not in exclude_types:
            imaging_types = IMAGING_TYPES.copy()
            imaging_types.update(['mir_lrs-fixedslit', 'mir_lrs-slitless'])
            if output_model.meta.exposure.type.lower() in imaging_types:
                try:
                    update_s_region_imaging(output_model)
                except Exception as exc:
                    log.error("Unable to update S_REGION for type {}: {}".format(
                        output_model.meta.exposure.type, exc))
                else:
                    log.info("assign_wcs updated S_REGION to {0}".format(
                        output_model.meta.wcsinfo.s_region))
                if output_model.meta.exposure.type.lower() == 'mir_lrs-slitless':
                    output_model.wavelength = get_wavelengths(output_model)
            elif output_model.meta.exposure.type.lower() == "nrs_ifu":
                update_s_region_nrs_ifu(output_model, mod)
            elif output_model.meta.exposure.type.lower() == 'mir_mrs':
                update_s_region_mrs(output_model)
            else:
                try:
                    update_s_region_spectral(output_model)
                except Exception as exc:
                    log.info("Unable to update S_REGION for type {}: {}".format(
                        output_model.meta.exposure.type, exc))

    log.info("COMPLETED assign_wcs")
    return output_model
