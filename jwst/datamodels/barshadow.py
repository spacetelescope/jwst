from .reference import ReferenceFileModel

__all__ = ['BarshadowModel']


class BarshadowModel(ReferenceFileModel):
    """
    A data model for Bar Shadow correction information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        Array defining the bar shadow correction as a function of Y and
        wavelength.

    variance : numpy array
        Variance array.

    """
    schema_url = "barshadow.schema.yaml"
