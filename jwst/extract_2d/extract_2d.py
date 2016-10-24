#
#  Module for 2d extraction
#
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging
import copy
import numpy as np
from astropy.modeling import models as astmodels
from .. import datamodels
from asdf import AsdfFile
from ..assign_wcs import nirspec


log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def extract2d(input_model, which_subarray=None):
    exp_type = input_model.meta.exposure.type.upper()
    log.info('EXP_TYPE is {0}'.format(exp_type))
    if exp_type in ['NRS_FIXEDSLIT', 'NRS_MSASPEC']:
        if which_subarray is None:
            open_slits = input_model.meta.wcs.get_transform('gwa', 'slit_frame').slits
        else:
            #open_slits = [nirspec.slit_name2id[which_subarray]]
            raise NotImplementedError()
    else:
        # Set the step status to SKIPPED, since it won't be done.
        input_model.meta.cal_step.extract_2d = 'SKIPPED'
        return input_model

    log.debug('open slits {0}'.format(open_slits))

    output_model = datamodels.MultiSlitModel()
    output_model.update(input_model)

    for slit in open_slits:
        slit_wcs = nirspec.nrs_wcs_set_input(input_model, slit.name)
        xlo, xhi = slit_wcs.domain[0]['lower'], slit_wcs.domain[0]['upper']
        ylo, yhi = slit_wcs.domain[1]['lower'], slit_wcs.domain[1]['upper']

        log.info('Name of subarray extracted: %s', slit.name)
        log.info('Subarray x-extents are: %s %s', xlo, xhi)
        log.info('Subarray y-extents are: %s %s', ylo, yhi)

        ext_data = input_model.data[ylo: yhi, xlo: xhi].copy()
        ext_err = input_model.err[ylo: yhi, xlo: xhi].copy()
        ext_dq = input_model.dq[ylo: yhi, xlo: xhi].copy()
        new_model = datamodels.ImageModel(data=ext_data, err=ext_err, dq=ext_dq)
        shape = ext_data.shape
        domain = [{'lower': 0, 'upper': shape[1], 'includes_lower': True, 'includes_upper': False},
                  {'lower':0, 'upper': shape[0], 'includes_lower': True, 'includes_upper': False}]
        slit_wcs.domain = domain
        new_model.meta.wcs = slit_wcs
        output_model.slits.append(new_model)
        # set x/ystart values relative to the image (screen) frame.
        # The overall subarray offset is recorded in model.meta.subarray.
        nslit = len(output_model.slits) - 1
        output_model.slits[nslit].name = str(slit.name)
        #output_model.slits[nslit].source_id = slit.source_id
        output_model.slits[nslit].xstart = xlo + 1
        output_model.slits[nslit].xsize = xhi - xlo
        output_model.slits[nslit].ystart = ylo + 1
        output_model.slits[nslit].ysize = yhi - ylo
    del input_model
    # Set the step status to COMPLETE
    output_model.meta.cal_step.extract_2d = 'COMPLETE'
    return output_model

