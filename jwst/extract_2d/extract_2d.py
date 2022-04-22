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
              reference_files={},
              grism_objects=None,
              tsgrism_extract_height=None,
              wfss_extract_half_height=None,
              extract_orders=None,
              mmag_extract=None,
              nbright=None):
    """
    The main extract_2d function

    Parameters
    ----------
    input_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodels.CubeModel`
        Input data model.
    slit_name : str or int
        Slit name.
    reference_files : dict
        Reference files.
    grism_objects : list
        A list of grism objects.
    tsgrism_extract_height: int
        Cross-dispersion extraction height to use for time series grisms.
        This will override the default which for NRC_TSGRISM is a set
        size of 64 pixels.
    wfss_extract_half_height : int
        Cross-dispersion extraction half height in pixels, WFSS mode.
        Overwrites the computed extraction height.
    extract_orders : list
        A list of spectral orders to be extracted.
    mmag_extract : float
        Minimum (faintest) abmag to extract for WFSS mode.
    nbright : float
        Number of brightest objects to extract, WFSS mode.

    Returns
    -------
    output_model : `~jwst.datamodels.ImageModel` or `~jwst.datamodelsCubeModel`
      A copy of the input_model that has been processed.

    """
    nrs_modes = ['NRS_FIXEDSLIT', 'NRS_MSASPEC', 'NRS_BRIGHTOBJ',
                 'NRS_LAMP', 'NRS_AUTOFLAT', 'NRS_AUTOWAVE']
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
            if tsgrism_extract_height is None:
                tsgrism_extract_height = 64
            output_model = extract_tso_object(input_model,
                                              reference_files=reference_files,
                                              tsgrism_extract_height=tsgrism_extract_height,
                                              extract_orders=extract_orders)
        else:
            output_model = extract_grism_objects(input_model,
                                                 grism_objects=grism_objects,
                                                 reference_files=reference_files,
                                                 extract_orders=extract_orders,
                                                 mmag_extract=mmag_extract,
                                                 wfss_extract_half_height=wfss_extract_half_height,
                                                 nbright=nbright)

    else:
        log.info(f'EXP_TYPE {exp_type} not supported for extract 2D')
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model

    # Set the step status to COMPLETE
    output_model.meta.cal_step.extract_2d = 'COMPLETE'
    del input_model
    return output_model
