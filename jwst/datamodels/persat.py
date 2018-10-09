from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['PersistenceSatModel']

class PersistenceSatModel(ReferenceFileModel):
    """
    A data model for the persistence saturation value (full well).

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
    schema_url = "persat.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(PersistenceSatModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
