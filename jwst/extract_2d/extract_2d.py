#
#  Top level module for 2d extraction.
#

import logging

from .nirspec import nrs_extract2d
from .grisms import extract_grism_objects, extract_tso_object

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def extract2d(input_model,
              slit_name=None,
              reference_file=None,
              grism_objects=None,
              extract_height=None,
              extract_orders=None,
              mmag_extract=99.):
    """
    The main extract_2d function

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
    slit_name : str or int
        Slit name.
    reference_file : str
        Reference file name.
    grism_objects : list
        A list of grism objects.
    extract_height: int
        Cross-dispersion extraction height to use for time series grisms.
        This will override the default which for NRC_TSGRISM is a set
        size of 64 pixels.

    Returns
    -------
    output_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodelsCubeModel`
      A copy of the input_model that has been processed.

    """
    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                 'NRS_LAMP', 'NRS_AUTOFLAT']
    slitless_modes = ['NIS_WFSS', 'NRC_WFSS', 'NRC_TSGRISM']

    exp_type = input_model.meta.exposure.type.upper()
    log.info(f'EXP_TYPE is {exp_type}')

    if exp_type in nrs_modes:
        if input_model.meta.instrument.grating.lower() == "mirror":
            # Catch the case of EXP_TYPE=NRS_LAMP and grating=MIRROR
            log.info(f'EXP_TYPE {exp_type} with grating=MIRROR not supported for extract 2D')
            input_model.meta.cal_step.extract_2d = 'SKIPPED'
            return input_model
        output_model = nrs_extract2d(input_model, slit_name=slit_name)
    elif exp_type in slitless_modes:
        if exp_type == 'NRC_TSGRISM':
            if extract_height is None:
                extract_height = 64
            output_model = extract_tso_object(input_model,
                                              reference_file=reference_file,
                                              extract_height=extract_height,
                                              extract_orders=extract_orders)
        else:
            output_model = extract_grism_objects(input_model,
                                                 grism_objects=grism_objects,
                                                 reference_file=reference_file,
                                                 extract_orders=extract_orders,
                                                 mmag_extract=99.)

    else:
        log.info(f'EXP_TYPE {exp_type} not supported for extract 2D')
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model

    # Set the step status to COMPLETE
    output_model.meta.cal_step.extract_2d = 'COMPLETE'
    del input_model
    return output_model
