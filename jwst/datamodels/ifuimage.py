from .image import ImageModel
from .model_base import DataModel


__all__ = ['IFUImageModel']


class IFUImageModel(DataModel):
    """
    A data model for 2D IFU images.

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

    relsens2d: numpy array
        The relative sensitivity 2D array.
    """
    schema_url = "ifuimage.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, ImageModel):
            super(IFUImageModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            if init.hasattr('area'):
                self.area = init.area
            if init.hasattr('relsens2d'):
                self.relsens2d = init.relsens2d
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            return

        super(IFUImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
