from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['FlatModel']


class FlatModel(ReferenceFileModel):
    """
    A data model for 2D flat-field images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  2-D.

    dq : numpy array
        The data quality array.  2-D.

    err : numpy array
        The error array.  2-D.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "flat.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(FlatModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
