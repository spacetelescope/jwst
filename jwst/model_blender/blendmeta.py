"""
Merge metadata from multiple models.

This module will create a new metadata instance and table from a
list of input datamodels or filenames.
"""

from .blender import ModelBlender


__all__ = ["blendmodels"]

# Primary functional interface for the code


def blendmodels(product, inputs, ignore=None):
    """
    Blend datamodel metadata.

    Parameters
    ----------
    product : JwstDataModel
        A datamodel that will have its metadata set
        to the blended metadata and have the metadata
        table assigned to the "hdrtab" attribute.

    inputs : list of JwstDataModel
        Input datamodels with metadata to blend.

    ignore : list of str, optional
        A list of metadata attributes to ignore during blending.
        These attributes will not be set on the output/combined.
        These attributes must be strings containing the dotted
        path of each attribute (for example "meta.filename").
        (Note that "meta.wcs" will always be ignored).
    """
    blender = ModelBlender(blend_ignore_attrs=ignore)
    for model in inputs:
        blender.accumulate(model)
    blender.finalize_model(product)
