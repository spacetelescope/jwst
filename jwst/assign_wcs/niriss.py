from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u

import gwcs.coordinate_frames as cf
from . import pointing
from .util import not_implemented_mode


def create_pipeline(input_model, reference_files):
    '''
    get reference files from crds

    '''
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)

    return pipeline


def imaging(input_model, reference_files):
    """
    The NIRISS imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'jwst_niriss_distortioon_0001.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    focal = cf.Frame2D(name='focal', axes_order=(0, 1), unit=(u.arcmin, u.arcmin))
    sky = cf.CelestialFrame(reference_frame=coord.ICRS())
    distortion = imaging_distortion(input_model, reference_files)
    fitswcs_transform = pointing.create_fitswcs_transform(input_model)
    pipeline = [(detector, distortion),
                (focal, None)]
                #(sky, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    return distortion


exp_type2transform = {'nis_image': imaging,
                      'nis_wfss': not_implemented_mode,#        ?? WFSS spec
                      'nis_soss': not_implemented_mode,#       ?? FS spec
                      'nis_ami': not_implemented_mode#        ?? imaging
                      }
