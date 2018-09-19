from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['FringeModel']


class FringeModel(ReferenceFileModel):
    """
    A data model for 2D fringe correction images.

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
    schema_url = "fringe.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(FringeModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
