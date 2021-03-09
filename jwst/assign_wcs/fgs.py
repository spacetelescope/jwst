"""
FGS WCS pipeline - depends on EXP_TYPE.
"""
import logging

from astropy import units as u
from astropy import coordinates as coord
from gwcs import coordinate_frames as cf

from .util import (not_implemented_mode, subarray_transform,
                   transform_bbox_from_shape, bounding_box_from_subarray)
from . import pointing
from ..datamodels import DistortionModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging"]


def create_pipeline(input_model, reference_files):
    """
    Create a ``gWCS.pipeline`` using models from reference files.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        The data model.
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

    It includes 3 coordinate frames - "detector", "v2v3" and "world".
    Uses a ``distortion`` reference file.

    Parameters
    ----------
    input_model : `~jwst.datamodels.DataModel`
        The data model.
    reference_files : dict
        {reftype: file_name} mapping.
        Reference files.

    Returns
    -------
    pipeline : list
        The WCS pipeline.
    """
    # Create coordinate frames for the ``imaging`` mode.
    detector = cf.Frame2D(name='detector', axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(name='v2v3', axes_order=(0, 1), axes_names=('v2', 'v3'),
                      unit=(u.arcsec, u.arcsec))
    v2v3vacorr = cf.Frame2D(name='v2v3vacorr', axes_order=(0, 1),
                            axes_names=('v2', 'v3'), unit=(u.arcsec, u.arcsec))
    world = cf.CelestialFrame(name='world', reference_frame=coord.ICRS())

    # Create the v2v3 to sky transform.
    tel2sky = pointing.v23tosky(input_model)

    # Read the distortion from the reference file
    distortion = imaging_distortion(input_model, reference_files)

    # If subarray, create an offset transform to be prepended to the distortion.
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        # Assign a bounding_box based on subarray's xsize and ysize
        distortion = subarray2full | distortion
        distortion.bounding_box = bounding_box_from_subarray(input_model)

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref
    )

    pipeline = [(detector, distortion),
                (v2v3, va_corr),
                (v2v3vacorr, tel2sky),
                (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the transform from "detector" to "v2v3".
    """
    dist = DistortionModel(reference_files['distortion'])
    transform = dist.model

    # Check if the transform in the reference file has a ``bounding_box``.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        transform.bounding_box
    except NotImplementedError:
        transform.bounding_box = transform_bbox_from_shape(input_model.data.shape)
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
