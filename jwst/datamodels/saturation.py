from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['SaturationModel']

class SaturationModel(ReferenceFileModel):
    """
    A data model for saturation checking information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "saturation.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(SaturationModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
