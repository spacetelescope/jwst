from .model_base import JwstDataModel


__all__ = ['SegmentationMapModel']


class SegmentationMapModel(JwstDataModel):
    """
    A data model for 2D segmentation maps

    Parameters
    __________
    data : numpy uint32 array
         The segmentation map
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/segmap.schema"
