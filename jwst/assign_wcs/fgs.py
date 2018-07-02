"""
FGS WCS pipeline - depends on EXP_TYPE.
"""
import logging

from astropy import units as u
from astropy import coordinates as coord
from gwcs import coordinate_frames as cf

from .util import not_implemented_mode, subarray_transform
from . import pointing
from ..datamodels import (DistortionModel, ImageModel, CubeModel)

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging"]


def create_pipeline(input_model, reference_files):
    """
    Create a ``gWCS.pipeline`` using models from reference files.

    Parameters
    ----------
    input_model : jwst.datamodels.DataModel
        Either an ImageModel or a CubeModel
    reference_files : dict
        {reftype: file_name} mapping.
        Reference files.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    log.info("Creating a FGS {0} pipeline with references {1}".format(
        exp_type, reference_files))
    return pipeline


def imaging(input_model, reference_files):
    """
    The FGS imaging WCS pipeline.

    It includes 3 coordinate frames -
    "detector", "v2v3" and "world".

    Uses a ``distortion`` reference file.
    """
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(name='world', reference_frame=coord.ICRS())
    # Crete the v2v3 to sky transform.
    tel2sky = pointing.v23tosky(input_model)

    # If subarray, ceate an offset transform to be prepended to the distortion.
    subarray2full = subarray_transform(input_model)
    if reference_files:
        imdistortion = imaging_distortion(input_model, reference_files)
        distortion = subarray2full | imdistortion
        # If the bounding box is saved in the model, move it to the first transform.
        distortion.bounding_box = imdistortion.bounding_box
        del imdistortion.bounding_box
    else:
        distortion = subarray2full

    pipeline = [(detector, distortion),
                (v2v3, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the transform from "detector" to "v2v3".
    """
    dist = DistortionModel(reference_files['distortion'])
    transform = dist.model

    # Get the ``bounding_box`` from the transform in the reference file.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        bb = transform.bounding_box
    except NotImplementedError:
        shape = input_model.data.shape
        # Note: Since bounding_box is attached to the model here
        # it's in reverse order.
        """
        A CubeModel is always treated as a stack (in dimension 1)
        of 2D images, as opposed to actual 3D data. In this case
        the bounding box is set to the 2nd and 3rd dimension.
        """
        if isinstance(input_model, CubeModel):
            bb = ((-0.5, shape[1] - 0.5),
                  (-0.5, shape[2] - 0.5))
        elif isinstance(input_model, ImageModel):
            bb = ((-0.5, shape[0] - 0.5),
                  (-0.5, shape[1] - 0.5))
        else:
            raise TypeError("Input is not an ImageModel or CubeModel")

        transform.bounding_box = bb
    dist.close()
    return transform


# EXP_TYPE to function mapping.
# The function creates the WCS pipeline.
exp_type2transform = {'fgs_image': imaging,
                      'fgs_focus': imaging,
                      'fgs_skyflat': not_implemented_mode,
                      'fgs_intflat': not_implemented_mode,
                      'fgs_dark': not_implemented_mode
                      }
