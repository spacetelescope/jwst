from .image import ImageModel
from .model_base import DataModel


__all__ = ['IFUImageModel']


class IFUImageModel(DataModel):
    """
    A data model for 2D IFU images.

    Attributes
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    area : numpy float32 array
         Pixel area map array

    relsens2d : numpy float32 array
         Sensitivity array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    wavelength : numpy float32 array
         wavelength

    pathloss_pointsource2d : numpy float32 array
         2-d array for pathloss (point source)

    pathloss_pointsource : numpy float32 array
         pathloss array for point sources

    wavelength_pointsource : numpy float32 array
         wavelength array for point sources

    pathloss_uniformsource2d : numpy float32 array
         2-d array for pathloss (uniform source)

    pathloss_uniformsource : numpy float32 array
         pathloss_array for uniform sources

    wavelength_uniformsource : numpy float32 array
         wavelength array for uniform sources
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
