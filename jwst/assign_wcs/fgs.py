"""
FGS WCS pipeline - depends on EXP_TYPE.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)
import logging

from asdf import AsdfFile
from astropy.modeling import models
from astropy import units as u
from astropy import coordinates as coord
from gwcs import coordinate_frames as cf

from .util import not_implemented_mode
from . import pointing


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def create_pipeline(input_model, reference_files):
    """
    Create a pipeline list based on EXP_TYPE.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Either an ImageModel or a CubeModel
    reference_files : dict
        {reftype: file_name} mapping
        In the pipeline it's returned by CRDS.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    log.info("Creating a FGS {0} pipeline with references {1}".format(
        exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    The FGS imaging pipeline includes 3 coordinate frames -
    detector, focal plane and sky.

    reference_files={'distortion': 'jwst_fgs_distortioon_0001.asdf'}
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(name='world', reference_frame=coord.ICRS())
    # V2, V3 to sky
    tel2sky = pointing.v23tosky(input_model)

    if reference_files:
        distortion = imaging_distortion(input_model, reference_files)
    else:
        distortion = models.Identity(2)
    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    distortion = AsdfFile.open(reference_files['distortion']).tree['model']
    # Convert to arcsec
    transform = distortion | models.Scale(60) & models.Scale(60)
    return transform


exp_type2transform = {'fgs_image': imaging,
                      'fgs_focus': imaging,
                      'fgs_skyflat': imaging,
                      'fgs_intflat': imaging,
                      'fgs_dark': not_implemented_mode
                      }
