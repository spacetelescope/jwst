from __future__ import absolute_import

import numpy as np
from astropy.io import fits
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

    # For NIRSpec MSA exposures, read the stellarity value for the
    # source in each extracted slit and set the point/extended value
    # based on the stellarity.
    elif exptype == 'NRS_MSASPEC':

        # Loop over the input slits
        for slit in input_model.slits:
            stellarity = slit.stellarity
            log.debug('source_id=%g, stellarity=%g' %
                      (slit.source_id, stellarity))

            # Eventually the stellarity value will be compared against
            # a threshold value from a reference file. For now, the
            # threshold is hardwired.
            if stellarity > 0.75:
                slit.source_type = 'POINT'
            else:
                slit.source_type = 'EXTENDED'

        # Set the source type value in the primary header to
        # a harmless default
        input_model.meta.target.source_type = 'UNKNOWN'

    # Unrecognized exposure type
    else:
        log.warning('EXP_TYPE %s not applicable to this operation' %
                    exptype)
        log.warning('Step will be skipped')
        return None

    # We're done
    return input_model
