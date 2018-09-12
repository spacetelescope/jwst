from .reference import ReferenceImageModel

__all__ = ['DarkMIRIModel']


class DarkMIRIModel(ReferenceImageModel):
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
