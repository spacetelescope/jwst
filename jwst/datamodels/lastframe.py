from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['LastFrameModel']


class LastFrameModel(ReferenceFileModel):
    """
    A data model for Last frame correction reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    err : numpy array
        The error array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "lastframe.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(LastFrameModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
