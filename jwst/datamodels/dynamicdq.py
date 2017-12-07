import numpy as np
from . import dqflags

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
log.addHandler(logging.NullHandler())

def dynamic_mask(input_model):
    #
    # Return a mask model given a mask with dynamic DQ flags
    # Dynamic flags define what each plane refers to using the DQ_DEF extension

    dq_table = input_model.dq_def
    # Get the DQ array and the flag definitions
    if (dq_table is not None and
        not np.isscalar(dq_table) and
        len(dq_table.shape) and
        len(dq_table)):
        #
        # Make an empty mask
        dqmask = np.zeros(input_model.dq.shape, dtype=input_model.dq.dtype)
        for record in dq_table:
            bitplane = record['VALUE']
            dqname = record['NAME'].strip()
            try:
                standard_bitvalue = dqflags.pixel[dqname]
            except KeyError:
                log.warning('Keyword %s does not correspond to an existing DQ mnemonic, so will be ignored' % (dqname))
                continue
            just_this_bit = np.bitwise_and(input_model.dq, bitplane)
            pixels = np.where(just_this_bit != 0)
            dqmask[pixels] = np.bitwise_or(dqmask[pixels], standard_bitvalue)
    else:
        dqmask = input_model.dq

    return dqmask
