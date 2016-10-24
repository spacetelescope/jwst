from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging
import importlib
from gwcs.wcs import WCS


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

    pipeline = mod.create_pipeline(input_model, reference_files)
    # Initialize the output model as a copy of the input
    # Make the copy after the WCS pipeline is created in order to pass updates to the model.
    output_model = input_model.copy()
    if pipeline is None:
        output_model.meta.cal_step.assign_wcs = 'SKIPPED'
    else:
        wcs = WCS(pipeline)
        #if output_model.meta.exposure.type.lower() in ['nrs_fixedslit', 'nrs_msaspec']:
            #record_open_slits(wcs)
        output_model.meta.wcs = wcs
        output_model.meta.cal_step.assign_wcs = 'COMPLETE'
    return output_model

