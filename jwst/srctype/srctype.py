from __future__ import absolute_import

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def set_source_type(input_model):
    """
    """

    # Get the exposure type of the input model
    try:
        exptype = input_model.meta.exposure.type
        log.debug('Input EXP_TYPE is %s' % exptype)
    except:
        log.error('Failed to access EXP_TYPE value in input')
        log.error('Step will be skipped')
        return None

    # For exposure types that use a single source, get the user-supplied
    # source type from the selection they provided in the APT
    if exptype in ['MIR_LRS-FIXEDSLIT', 'MIR_LRS-SLITLESS', 'MIR_MRS',
                   'NIS_SOSS', 'NRS_FIXEDSLIT', 'NRS_IFU']:

        # Get the value the user specified (if any)
        user_type = input_model.meta.target.source_type

        if (user_type is not None) and (user_type in ['POINT', 'EXTENDED']):

            # Use the value supplied by the user
            log.info('Using input SRCTYPE of %s' % user_type)
            src_type = user_type

        else:

            # Set a default value based on the exposure type
            if exptype == 'MIR_MRS':
                src_type = 'EXTENDED'
            else:
                src_type = 'POINT'
            log.info('Input SRCTYPE is unknown. Setting to default ' +
                     'value of %s' % src_type)

        input_model.meta.target.source_type = src_type
    
    return input_model
