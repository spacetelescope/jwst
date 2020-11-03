from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['DarkModel']


class DarkModel(ReferenceFileModel):
    """
    A data model for dark reference files.

    Parameters
    __________
    data : numpy float32 array
         Dark current array

    dq : numpy uint16 array
         2-D data quality array for all planes

    err : numpy float32 array
         Error array

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/dark.schema"

    def __init__(self, init=None, **kwargs):
        super(DarkModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
