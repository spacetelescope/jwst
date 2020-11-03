from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['SaturationModel']


class SaturationModel(ReferenceFileModel):
    """
    A data model for saturation checking information.

    Parameters
    __________
    data : numpy float32 array
         Saturation threshold

    dq : numpy uint32 array
         2-D data quality array for all planes

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/saturation.schema"

    def __init__(self, init=None, **kwargs):
        super(SaturationModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
