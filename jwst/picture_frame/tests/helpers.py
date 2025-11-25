import numpy as np
from stdatamodels.jwst import datamodels

from jwst.picture_frame import picture_frame as pf

__all__ = ["picture_frame_model"]


def picture_frame_model():
    """
    Create a picture frame datamodel.

    Returns
    -------
    pctfrm : `~stdatamodels.jwst.datamodels.PictureFrameModel`
        The picture frame reference file datamodel.
    """
    pctfrm = datamodels.PictureFrameModel()
    pctfrm.data = np.full((2048, 2048), 0.1)
    pctfrm.data[
        pf.CENTER_REGION[0] : pf.CENTER_REGION[1], pf.CENTER_REGION[0] : pf.CENTER_REGION[1]
    ] = 1.1

    # add required metadata
    pctfrm.meta.description = "test"
    pctfrm.meta.reftype = "test"
    pctfrm.meta.author = "test"
    pctfrm.meta.pedigree = "test"
    pctfrm.meta.useafter = "test"

    return pctfrm
