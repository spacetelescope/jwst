from stcal.dynamicdq import dynamic_mask
from .dqflags import pixel
from .reference import ReferenceFileModel


__all__ = ['LastFrameModel']


class LastFrameModel(ReferenceFileModel):
    """
    A data model for Last frame correction reference files.

    Parameters
    __________
    data : numpy float32 array
         Last Frame Correction array

    dq : numpy uint32 array
         2-D data quality array

    err : numpy float32 array
         Error array

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/lastframe.schema"

    def __init__(self, init=None, **kwargs):
        super(LastFrameModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self, pixel)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
