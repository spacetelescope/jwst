from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging
import importlib
from gwcs.wcs import WCS
from.util import update_s_region

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def load_wcs(input_model, reference_files={}):
    """
    Create the WCS object from the input model and reference files and
    store the pickled WCS-related object into the model meta data.
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

    # Add WCS keywords for the spectral axis.
    if input_model.meta.wcsinfo.wcsaxes == 3:
        _add_3rd_axis(input_model)

    pipeline = mod.create_pipeline(input_model, reference_files)

    # Initialize the output model as a copy of the input
    # Make the copy after the WCS pipeline is created in order to pass updates to the model.
    if pipeline is None:
        input_model.meta.cal_step.assign_wcs = 'SKIPPED'
        log.warning("SKIPPED assign_wcs")
        return input_model
    else:
        output_model = input_model.copy()
        wcs = WCS(pipeline)
        output_model.meta.wcs = wcs
        output_model.meta.cal_step.assign_wcs = 'COMPLETE'
        exclude_types = ['NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ', 'NRS_IFU',
                         'NRS_MSASPEC', 'NRS_LAMP', 'MIR_MRS', 'NIS_SOSS',
                         'NRC_GRISM', 'NRC_TSGRISM', 'NIS_WFSS']

        if output_model.meta.exposure.type not in exclude_types:
            orig_s_region = output_model.meta.wcsinfo.s_region.strip()
            update_s_region(output_model)
            if orig_s_region != output_model.meta.wcsinfo.s_region.strip():
                log.info("assign_wcs updated S_REGION to {0}".format(output_model.meta.wcsinfo.s_region))
        else:
            log.info("assign_wcs did not update S_REGION for type {0}".format(output_model.meta.exposure.type))
        log.info("COMPLETED assign_wcs")
    return output_model


def _add_3rd_axis(input_model):
    """
    Add WCS keywords and their default values for the spectral axis.

    SDP adds CTYPE3 and CUNIT3.

    """
    input_model.meta.wcsinfo.pc1_3 = 0.
    input_model.meta.wcsinfo.pc2_3 = 0.
    input_model.meta.wcsinfo.pc3_3 = 1.
    input_model.meta.wcsinfo.crval3 = 0.
    input_model.meta.wcsinfo.crpix3 = 0.
    input_model.meta.wcsinfo.cdelt3 = 1.
