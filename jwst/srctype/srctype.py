import logging
from ..lib import pipe_utils

from .. import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def set_source_type(input_model, source_type=None):
    """
    Set source_type based on APT input, user specification, exposure type,
    or default values.

    Parameters
    ----------
    input_model : `~jwst.datamodels.CubeModel`, `~jwst.datamodels.ImageModel`,
                  `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.MultiSlitModel`,
                  or `~jwst.datamodels.SlitModel`
        The data model to be processed.

    source_type : str,  {POINT, EXTENDED}
        User-requested value for source type.

    Returns
    -------
    input_model : `~jwst.datamodels.CubeModel`, `~jwst.datamodels.ImageModel`,
                  `~jwst.datamodels.IFUImageModel`, `~jwst.datamodels.MultiSlitModel`,
                  or `~jwst.datamodels.SlitModel`
        The updated model.
    """

    # Get the exposure type of the input model
    exptype = input_model.meta.exposure.type
    if exptype is None:
        log.error('EXP_TYPE value not found in input')
        raise RuntimeError('Step cannot be executed without an EXP_TYPE value')
    else:
        log.info(f'Input EXP_TYPE is {exptype}')
    # For exposure types that have a single source specification, get the
    # user-supplied source type from the selection they provided in the APT
    if exptype in ['MIR_LRS-FIXEDSLIT', 'MIR_LRS-SLITLESS', 'MIR_MRS',
                   'NRC_TSGRISM', 'NIS_SOSS', 'NRS_FIXEDSLIT',
                   'NRS_BRIGHTOBJ', 'NRS_IFU']:

        # Get info about the exposure, including whether it's a background
        # target and the dither pattern type
        bkg_target = input_model.meta.observation.bkgdtarg
        if exptype == 'MIR_MRS':
            patttype = input_model.meta.dither.optimized_for
        else:
            patttype = input_model.meta.dither.primary_type

        # The keyword SRCTYAPT was added in JWSTKD-354 to retain the value
        # given by the user in the APT, while the cal code then sets a value
        # for the SRCTYPE keyword for use in calibration. Try to preserve
        # backwards compatibility with old datasets by first looking for the
        # SRCTYAPT keyword and using it if available, and if not, then use
        # SRCTYPE as both input and output (as before).
        try:
            user_type = input_model.meta.target.source_type_apt
            log.info(f'Input SRCTYAPT = {user_type}')
            if user_type is None:
                log.warning('SRCTYAPT keyword not found in input; using SRCTYPE instead')
                user_type = input_model.meta.target.source_type
                input_model.meta.target.source_type_apt = user_type
        except AttributeError:
            log.warning('SRCTYAPT keyword not found in input; using SRCTYPE instead')
            user_type = input_model.meta.target.source_type
            input_model.meta.target.source_type_apt = user_type

        # Check to see if the user specified a source type
        if source_type is not None:
            source_type = str(source_type).upper()

            # Check if the exposure type is a mode that allows setting
            if exptype in ['MIR_LRS-FIXEDSLIT', 'MIR_LRS-SLITLESS',
                           'MIR_MRS', 'NRC_TSGRISM', 'NRS_FIXEDSLIT',
                           'NRS_BRIGHTOBJ', 'NRS_IFU']:

                src_type = source_type

            log.warning(f'Based on user-input, setting SRCTYPE = {src_type}')
            input_model.meta.target.source_type = src_type

        elif bkg_target:

            # If this image is flagged as a BACKGROUND target, set the
            # source type to EXTENDED regardless of any other settings
            src_type = 'EXTENDED'
            log.info(f'Exposure is a background target; setting SRCTYPE = {src_type}')

        elif pipe_utils.is_tso(input_model):

            # Treat all TSO exposures as a point source
            src_type = 'POINT'
            log.info(f'Input is a TSO exposure; setting SRCTYPE = {src_type}')

        elif (patttype is not None) and (('NOD' in patttype) or ('POINT-SOURCE' in patttype)):

            # Set all nodded exposures to POINT source type
            src_type = 'POINT'
            log.info(f'Exposure is nodded; setting SRCTYPE = {src_type}')

        elif user_type in ['POINT', 'EXTENDED']:

            # Use the value supplied by the user
            src_type = user_type
            log.info(f'Using input source type = {src_type}')

        else:

            # Set a default value based on the exposure type
            if exptype in ['MIR_MRS', 'NRS_IFU']:
                src_type = 'EXTENDED'
            else:
                src_type = 'POINT'

            log.info(f'Input source type is unknown; setting default SRCTYPE = {src_type}')

        # Set the source type in the global meta attribute
        input_model.meta.target.source_type = src_type

        # If the input contains one or more slit instances,
        # set the value in each slit too
        if isinstance(input_model, datamodels.SlitModel):
            input_model.source_type = src_type

        elif input_model.meta.exposure.type == 'NRS_FIXEDSLIT':

            # NIRSpec fixed-slit is a special case: Apply the source type
            # determined above to only the primary slit (the one in which
            # the target is located). Set all other slits to the default
            # value, which for NRS_FIXEDSLIT is 'POINT'.
            default_type = 'EXTENDED'
            primary_slit = input_model.meta.instrument.fixed_slit
            log.debug(f' primary_slit = {primary_slit}')
            for slit in input_model.slits:
                if slit.name == primary_slit:
                    slit.source_type = src_type
                else:
                    slit.source_type = default_type
                log.debug(f' slit {slit.name} = {slit.source_type}')

    # For NIRSpec MSA exposures, read the stellarity value for the
    # source in each extracted slit and set the point/extended value
    # based on the stellarity.
    elif exptype == 'NRS_MSASPEC':

        # Loop over the input slits
        for slit in input_model.slits:
            stellarity = slit.stellarity

            # Eventually the stellarity value will be compared against
            # a threshold value from a reference file. For now, the
            # threshold is hardwired.
            if stellarity < 0.0:
                slit.source_type = 'POINT'
            elif stellarity > 0.75:
                slit.source_type = 'POINT'
            else:
                slit.source_type = 'EXTENDED'

            log.info(f'source_id={slit.source_id}, stellarity={stellarity:.4f}, type={slit.source_type}')

        # Remove the global target source type, so that it never mistakenly
        # gets used for MOS data, which should always use slit-specific values
        input_model.meta.target.source_type = None

    # Set all TSO exposures to POINT
    elif pipe_utils.is_tso(input_model):
        src_type = 'POINT'
        log.info(f'Input is a TSO exposure; setting default SRCTYPE = {src_type}')
        input_model.meta.target.source_type = src_type

    # For WFSS modes check slit values of is_extended to set SRCTYPE
    elif exptype in ['NIS_WFSS', 'NRC_WFSS']:
        for slit in input_model.slits:
            if slit.is_extended:
                slit.source_type = 'EXTENDED'
            else:
                slit.source_type = 'POINT'
            log.info(f'source_id={slit.source_id}, type={slit.source_type}')

    # Unrecognized exposure type; set to UNKNOWN as default
    else:
        log.warning(f'EXP_TYPE {exptype} not applicable to this operation')
        src_type = 'UNKNOWN'
        log.warning(f'Setting SRCTYPE = {src_type}')
        input_model.meta.target.source_type = src_type

    # We're done
    return input_model
