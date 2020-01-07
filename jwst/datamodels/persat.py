from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['PersistenceSatModel']

class PersistenceSatModel(ReferenceFileModel):
    """
    A data model for the persistence saturation value (full well).

    Parameters
    __________
    data : numpy float32 array
         Persistence saturation threshold

    dq : numpy uint32 array
         data quality array

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/persat.schema"

    def __init__(self, init=None, **kwargs):
        super(PersistenceSatModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
