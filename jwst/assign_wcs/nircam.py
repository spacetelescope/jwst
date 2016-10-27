from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from asdf import AsdfFile
from astropy import coordinates as coord
from astropy import units as u
from astropy.modeling.models import Scale

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
    The NIRCAM imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky

    reference_files={'distortion': 'test.asdf', 'filter_offsets': 'filter_offsets.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.deg, u.deg))
    world = cf.CelestialFrame(reference_frame=coord.ICRS(), name='world')
    distortion = imaging_distortion(input_model, reference_files)
    tel2sky = pointing.v23tosky(input_model)
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    # Convert to arcsec
    transform = distortion | Scale(1/60) & Scale(1/60)
    return transform


exp_type2transform = {'nrc_image': imaging,
                      'nrc_slitless': not_implemented_mode, #WFSS mode
                      'nrc_tacq': not_implemented_mode,#       ?? distortion
                      'nrc_coron': not_implemented_mode,#    ?? distortion
                      'nrc_focus': not_implemented_mode,#       ?? distortion
                      'nrc_tss': not_implemented_mode,# custom SOSS like mode TBC
                      'nrc_tsi': not_implemented_mode,# custom soss like mode (TBC
                      'nrc_led': not_implemented_mode,#    ?? WFSS mode
                      }
