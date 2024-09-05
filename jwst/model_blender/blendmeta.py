"""blendmeta - Merge metadata from multiple models.

    This module will create a new metadata instance and table from a
    list of input datamodels or filenames.
"""

from .blender import ModelBlender


# Primary functional interface for the code

def blendmodels(product, inputs, ignore=None):
    blender = ModelBlender(blend_ignore_attrs=ignore)
    for model in inputs:
        blender.accumulate(model)
    blender.finalize_model(product)
