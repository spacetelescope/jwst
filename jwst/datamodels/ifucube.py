from .model_base import DataModel


__all__ = ['IFUCubeModel']


class IFUCubeModel(DataModel):
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

    def __init__(self, init=None, **kwargs):
        super(IFUCubeModel, self).__init__(init=init, **kwargs)

       # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
