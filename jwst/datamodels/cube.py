from . import model_base


__all__ = ['CubeModel']


class CubeModel(model_base.DataModel):
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

    def __init__(self, init=None, data=None, dq=None, err=None, zeroframe=None,
                 relsens=None, int_times=None, area=None, wavelength=None,
                 var_poisson=None, var_rnoise=None, **kwargs):
      
        super(CubeModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if relsens is not None:
            self.relsens = relsens

        if int_times is not None:
            self.int_times = int_times

        if area is not None:
            self.area = area

        if wavelength is not None:
            self.wavelength = wavelength

        if var_poisson is not None:
            self.var_poisson = var_poisson

        if var_rnoise is not None:
            self.var_rnoise = var_rnoise
   
        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
