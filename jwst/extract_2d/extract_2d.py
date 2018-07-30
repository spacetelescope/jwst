#
#  Top level module for 2d extraction.
#

import logging

from .nirspec import nrs_extract2d
from .grisms import extract_grism_objects

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract2d(input_model, slit_name=None, apply_wavecorr=False, reference_files={}, grism_objects=[]):
    """
    The main extract_2d function.

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
    slit_name : str or int
        Slit name.
    apply_wavecorr : bool
        Flag whether to apply the zero point wavelength correction to
        Nirspec exposures.
    reference_files : dict
        Reference files.
    grism_objects : list
        A list of grism objects.
    """

    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                 'NRS_LAMP', 'NRS_AUTOFLAT']
    grism_modes = ['NIS_WFSS', 'NRC_WFSS']

    exp_type = input_model.meta.exposure.type.upper()
    log.info('EXP_TYPE is {0}'.format(exp_type))

    if exp_type in nrs_modes:
        output_model = nrs_extract2d(input_model, slit_name=slit_name,
                                     apply_wavecorr=apply_wavecorr,
                                     reference_files=reference_files)

    elif exp_type in grism_modes:
        output_model = extract_grism_objects(input_model, grism_objects=grism_objects, reference_files=reference_files)

    else:
        log.info("'EXP_TYPE {} not supported for extract 2D".format(exp_type))
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model


    # Set the step status to COMPLETE
    output_model.meta.cal_step.extract_2d = 'COMPLETE'
    del input_model
    return output_model
