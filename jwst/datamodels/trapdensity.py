from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['TrapDensityModel']

class TrapDensityModel(ReferenceFileModel):
    """
    A data model for the trap density of a detector, for persistence.

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data.

    dq: numpy array
        The data quality array.

    dq_def: numpy array
        The data quality definitions table.
    """
    schema_url = "trapdensity.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, dq_def=None, **kwargs):
        super(TrapDensityModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if dq_def is not None:
            self.dq_def = dq_def

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
