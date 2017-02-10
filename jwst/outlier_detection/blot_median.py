"""
Create a blotted images that match the inputs from the combined median image.

:Authors: Warren Hack

:License:

"""
import numpy as np

from ..resample import gwcs_blot
from .. import datamodels


def do_blot(median_model, input_models, **pars):
    # start by interpreting input parameters
    interp = pars.get('interp', 'poly5')
    sinscl = pars.get('sinscl', 1.0)

    # Initialize output product
    blot_models = datamodels.ModelContainer()

    blot = gwcs_blot.GWCSBlot(median_model)

    for model in input_models:
        blotted_median = model.copy()
        blot_root = '_'.join(model.meta.filename.replace('.fits', '').split('_')[:-1])
        blotted_median.meta.filename = '{}_blot.fits'.format(blot_root)

        # clean out extra data not related to blot result
        blotted_median.err = None
        blotted_median.dq = None
        # apply blot to re-create model.data from median image
        blotted_median.data = blot.extract_image(model, interp=interp,
            sinscl=sinscl)
        blot_models.append(blotted_median)
    return blot_models
