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
    interp = pars.get('interp','poly5')
    sinscl = pars.get('sinscl',1.0)
    
    # Initialize output product
    blot_models = datamodels.ModelContainer()

    do_blot = gwcs_blot.GWCSBlot(median_model)
            
    for input_img in input_models:
        blot_model = input_img.copy()
        blot_root = '_'.join(input_img.meta.filename.replace('.fits','').split('_')[:-1])
        blot_model.meta.filename = '{}_blot.fits'.format(blot_root)

        # clean out extra data not related to blot result
        blot_model.err = None
        blot_model.dq = None
        # apply blot to re-create input_img.data from median image
        blot_model.data = do_blot.extract_image(input_img.meta.wcs,
                                        interp=interp,sinscl=sinscl)
        blot_models.append(blot_model)
    return blot_models
    
