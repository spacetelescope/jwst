from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['DarkMIRIModel']


class DarkMIRIModel(ReferenceFileModel):
    """
    A data model for dark MIRI reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data (integration dependent)

    dq : numpy array
        The data quality array. (integration dependent)

    err : numpy array (integration dependent)
        The error array.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "darkMIRI.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(DarkMIRIModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
