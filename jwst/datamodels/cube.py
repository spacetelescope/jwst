from .model_base import DataModel


__all__ = ['CubeModel']


class CubeModel(DataModel):
    """
    A data model for 3D image cubes.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.  3-D.

    dq : numpy array
        The data quality array.  3-D.

    err : numpy array
        The error array.  3-D

    zeroframe: numpy array
        The zero-frame array.  3-D

    relsens: numpy array
        The relative sensitivity array.

    int_times : table
        The int_times table

    area: numpy array
        The pixel area array.  2-D

    wavelength: numpy array
        The wavelength array.  2-D

    var_poisson: numpy array
        The variance due to Poisson noise array.  3-D

    var_rnoise: numpy array
        The variance due to read noise array.  3-D
    """
    schema_url = "cube.schema.yaml"

    def __init__(self, init=None, **kwargs):

        super(CubeModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

