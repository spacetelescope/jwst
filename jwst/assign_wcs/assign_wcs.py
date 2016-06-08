from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import logging
import os
import importlib

import numpy as np
from astropy.io import fits
from gwcs.wcs import WCS

from .util import is_fits
from . import pointing


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def load_wcs (input_model, reference_files={}):
    """
    Create the WCS object from the input model and reference files and
    store the pickled WCS-related object into the model meta data.
    """
    exp_type = input_model.meta.exposure.type.lower()
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
    mod = importlib.import_module('.' + instrument, 'jwst_pipeline.assign_wcs')

    pipeline = mod.create_pipeline(input_model, reference_files)
    # Initialize the output model as a copy of the input
    # Make the copy after the WCS pipeline is created in order to pass updates to the model.
    output_model = input_model.copy()
    if pipeline is None:
        output_model.meta.cal_step.assign_wcs = 'SKIPPED'
    else:
        wcs = WCS(pipeline)
        output_model.meta.wcs = wcs
        output_model.meta.cal_step.assign_wcs = 'COMPLETE'
    return output_model


# extract2d
'''
def nrs_set_inputs(input_model):
    from . import nirspec
    from asdf import AsdfFile

    # get the file name from the model
    msa_status = "SPCB-GD-A.msa.fits"
    open_slits_id = nirspec.get_open_msa_slits(msa_status)

    filter = input_model.meta.instrument.filter
    grating = input_model.meta.instrument.grating
    frange = input_model.meta.ref_file.wavelengthrange.name
    wave_range = AsdfFile.open(frange)

    try:
        wrange = wave_range.tree['filter_grating'][filter + '_' + grating]['range']
    except KeyError:
        log.warning("Combination {0} missing in wavelengthrange file, cannot computs slit domain.".format(
            filter+'_'+grating))
    #for slit in open_slits_id:
    slit_wcs = nirspec.nrsmsa_wcs_set_input(input_model.meta.wcs, slit[0], slit[1], wrange)

    # input_model here is the MultiSlit model
    # attach the slit_wcs object to every slit.
    # at the end delete the main wcs object.
    wave_range.close()
    return slit_wcs


def nrs_ifu_set_inputs(input_model):
    from . import nirspec
    from asdf import AsdfFile
    slits = np.arange(30)

    filter = input_model.meta.instrument.filter
    grating = input_model.meta.instrument.grating
    frange = input_model.meta.ref_file.wavelengthrange.name
    wave_range = AsdfFile.open(frange)

    try:
        wrange = wave_range.tree['filter_grating'][filter + '_' + grating]['range']
    except KeyError:
        log.warning("Combination {0} missing in wavelengthrange file, cannot computs slit domain.".format(
            filter+'_'+grating))
    for slit in slits:
        slit_wcs = nirspec.nrs_wcs_set_input(input_model.meta.wcs, 0, slit, wrange)
    # input_model here is the MultiSlit model
    # attach the slit_wcs object to every slit.
    # at the end delete the main wcs object. 
    wave_range.close()
    return slit_wcs
'''
