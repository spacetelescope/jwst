from .model_base import JwstDataModel
from .image import ImageModel


__all__ = ['SlitModel', 'SlitDataModel']


class SlitDataModel(JwstDataModel):
    """
    A data model for 2D slit images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    wavelength : numpy float32 array
         Wavelength array, corrected for zero-point

    barshadow : numpy float32 array
         Bar shadow correction

    flatfield_point : numpy float32 array
         flatfield array for point source

    flatfield_uniform : numpy float32 array
         flatfield array for uniform source

    pathloss_point : numpy float32 array
         pathloss array for point source

    pathloss_uniform : numpy float32 array
         pathloss array for uniform source

    photom_point : numpy float32 array
         photom array for point source

    photom_uniform : numpy float32 array
         photom array for uniform source

    area : numpy float32 array
         Pixel area map array
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/slitdata.schema"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitDataModel, self).__init__(init=None, **kwargs)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.area = init.area
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
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


class SlitModel(JwstDataModel):
    """
    A data model for 2D images.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    wavelength : numpy float32 array
         Wavelength array, corrected for zero-point

    barshadow : numpy float32 array
         Bar shadow correction

    flatfield_point : numpy float32 array
         flatfield array for point source

    flatfield_uniform : numpy float32 array
         flatfield array for uniform source

    pathloss_point : numpy float32 array
         pathloss array for point source

    pathloss_uniform : numpy float32 array
         pathloss array for uniform source

    photom_point : numpy float32 array
         photom array for point source

    photom_uniform : numpy float32 array
         photom array for uniform source

    area : numpy float32 array
         Pixel area map array

    int_times : numpy table
         table of times for each integration
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/slit.schema"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(SlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.data = init.data
            self.dq = init.dq
            self.err = init.err
            self.area = init.area
            if init.hasattr('var_poisson'):
                self.var_poisson = init.var_poisson
            if init.hasattr('var_rnoise'):
                self.var_rnoise = init.var_rnoise
            if init.hasattr('wavelength'):
                self.wavelength = init.wavelength
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
