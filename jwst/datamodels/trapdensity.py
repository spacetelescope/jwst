from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['TrapDensityModel']

class TrapDensityModel(ReferenceFileModel):
    """
    A data model for the trap density of a detector, for persistence.

    Parameters
    __________
    data : numpy float32 array
         Trap density

    dq : numpy uint32 array
         data quality array

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/trapdensity.schema"

    def __init__(self, init=None, **kwargs):
        super(TrapDensityModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
