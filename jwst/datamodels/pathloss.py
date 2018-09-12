from .reference import ReferenceFileModel


__all__ = ['PathlossModel']


class PathlossModel(ReferenceFileModel):
    """
    A data model for pathloss correction information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    pointsource : numpy array
        Array defining the pathloss parameter for point sources.

    psvar : numpy array
        Variance array.

    uniform : numpy array
        Pathloss parameter for uniform illumination
    """
    schema_url = "pathloss.schema.yaml"
