from .model_base import DataModel
from .image import ImageModel


__all__ = ['SlitModel', 'SlitDataModel']


class SlitDataModel(DataModel):
    """
    A data model for 2D images.

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

    relsens : numpy array
        The relative sensitivity table.

    """

    schema_url = "slitdata.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitDataModel, self).__init__(init=None, **kwargs)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            for key in kwargs:
                setattr(self, key, kwargs[key])

            if init.meta.hasattr('wcs'):
                self.meta.wcs = init.meta.wcs
            else:
                self.meta.wcs = None
        else:
            super(SlitDataModel, self).__init__(init=init, **kwargs)
            if kwargs:
                for key in kwargs:
                    setattr(self, key, kwargs[key])


class SlitModel(DataModel):
    """
    A data model for 2D images.

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

    relsens : numpy array
        The relative sensitivity table.

    int_times : table
        The int_times table

    """
    schema_url = "slit.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.relsens = init.relsens
            self.area = init.area
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            if init.hasattr('int_times'):
                self.int_times = init.int_times
            if init.meta.hasattr('wcs'):
                self.meta.wcs = init.meta.wcs
            else:
                self.meta.wcs = None
        else:
            super(SlitModel, self).__init__(init=init, **kwargs)
            if kwargs:
                for key in kwargs:
                    setattr(self, key, kwargs[key])
