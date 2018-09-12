from .image import ImageModel


__all__ = ['IFUCubeModel']


class IFUCubeModel(ImageModel):
    """
    A data model for 3D IFU  cubes.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data.  3-D.

    dq: numpy array
        The data quality array.  3-D.

    err: numpy array
        The error array.  3-D

    weightmap: numpy array
        The weight map array.  3-D

    wavetable:  1-D table
        Optional table of  wavelengths of IFUCube slices

    """
    schema_url = "ifucube.schema.yaml"
