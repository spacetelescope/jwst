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

    # For NIRSpec MSA exposures, load the stellarity values from the
    # MSA configuration file that goes with this exposure
    if exptype == 'NRS_MSASPEC':

        # Get the MSA configuration file name and open the file
        msa_name = input_model.meta.instrument.msa_configuration_file
        msa_file = fits.open(msa_name)

        # Load the Source info table from the MSA file
        source_table = msa_file['SOURCE_INFO'].data

        # See if the stellarity column exists
        try:
            stellarities = source_table['stellarity']
        except:
            log.error('Stellarity data not found in MSA source table')
            log.error('Step will be skipped')
            msa_file.close()
            return None

        # Loop over the input slits and find matching source_id in table
        for slit in range(len(input_model.slits)):
            source_id = input_model.slits[slit].source_id
            row = np.where(source_table['source_id'] == source_id)

            # Get the stellarity value for this source
            stellarity = np.float32(source_table['stellarity'][row][0])
            log.debug('source_id=%g, stellarity=%g' % (source_id, stellarity))

            # Save the stellarity value in the slit meta data
            input_model.slits[slit].stellarity = stellarity

        # We're done; close the MSA configuration file
        msa_file.close()


    return input_model
