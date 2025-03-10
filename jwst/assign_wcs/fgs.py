"""FGS WCS pipeline - depends on EXP_TYPE."""

import logging

from astropy import units as u
from astropy import coordinates as coord
from astropy.modeling import bind_bounding_box
from gwcs import coordinate_frames as cf

from stdatamodels.jwst.datamodels import DistortionModel

from .util import (
    not_implemented_mode,
    subarray_transform,
    transform_bbox_from_shape,
    bounding_box_from_subarray,
)
from . import pointing

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


__all__ = ["create_pipeline", "imaging"]


def create_pipeline(input_model, reference_files):
    """
    Create a ``gWCS.pipeline`` using models from reference files.

    Parameters
    ----------
    input_model : `~jwst.datamodels.JwstDataModel`
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    exp_type = input_model.meta.exposure.type.lower()
    pipeline = exp_type2transform[exp_type](input_model, reference_files)
    log.info(f"Creating a FGS {exp_type} pipeline with references {reference_files}")
    return pipeline


def imaging(input_model, reference_files):
    """
    Create the WCS pipeline for FGS imaging data.

    It includes 3 coordinate frames - "detector", "v2v3" and "world".

    Parameters
    ----------
    input_model : ImageModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' reference files

    Returns
    -------
    pipeline : list
        The WCS pipeline, suitable for input into `gwcs.WCS`.
    """
    # Create coordinate frames for the ``imaging`` mode.
    detector = cf.Frame2D(name="detector", axes_order=(0, 1), unit=(u.pix, u.pix))
    v2v3 = cf.Frame2D(
        name="v2v3", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    v2v3vacorr = cf.Frame2D(
        name="v2v3vacorr", axes_order=(0, 1), axes_names=("v2", "v3"), unit=(u.arcsec, u.arcsec)
    )
    world = cf.CelestialFrame(name="world", reference_frame=coord.ICRS())

    # Create the v2v3 to sky transform.
    tel2sky = pointing.v23tosky(input_model)

    # Read the distortion from the reference file
    distortion = imaging_distortion(input_model, reference_files)

    # If subarray, create an offset transform to be prepended to the distortion.
    subarray2full = subarray_transform(input_model)
    if subarray2full is not None:
        # Assign a bounding_box based on subarray's xsize and ysize
        distortion = subarray2full | distortion

        # Bind the bounding box to the distortion model using the bounding box ordering
        # used by GWCS. This makes it clear the bounding box is set correctly to GWCS
        bind_bounding_box(distortion, bounding_box_from_subarray(input_model, order="F"), order="F")

    # Compute differential velocity aberration (DVA) correction:
    va_corr = pointing.dva_corr_model(
        va_scale=input_model.meta.velocity_aberration.scale_factor,
        v2_ref=input_model.meta.wcsinfo.v2_ref,
        v3_ref=input_model.meta.wcsinfo.v3_ref,
    )

    pipeline = [(detector, distortion), (v2v3, va_corr), (v2v3vacorr, tel2sky), (world, None)]
    return pipeline


def imaging_distortion(input_model, reference_files):
    """
    Create the transform from "detector" to "v2v3".

    Parameters
    ----------
    input_model : ImageModel
        The input data model.
    reference_files : dict
        Mapping between reftype (keys) and reference file name (vals).
        Requires 'distortion' reference file.

    Returns
    -------
    transform : `astropy.modeling.Model
        The transform from "detector" to "v2v3".
    """
    dist = DistortionModel(reference_files["distortion"])
    transform = dist.model

    # Check if the transform in the reference file has a ``bounding_box``.
    # If not set a ``bounding_box`` equal to the size of the image.
    try:
        transform.bounding_box  # noqa: B018
    except NotImplementedError:
        bind_bounding_box(
            transform, transform_bbox_from_shape(input_model.data.shape, order="F"), order="F"
        )

    dist.close()
    return transform


# EXP_TYPE to function mapping.
# The function creates the WCS pipeline.
exp_type2transform = {
    "fgs_image": imaging,
    "fgs_focus": imaging,
    "fgs_skyflat": not_implemented_mode,
    "fgs_intflat": not_implemented_mode,
    "fgs_dark": not_implemented_mode,
}
