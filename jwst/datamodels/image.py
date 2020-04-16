from .model_base import DataModel


__all__ = ['ImageModel']


class ImageModel(DataModel):
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

    zeroframe : numpy float32 array
         Zeroframe array

    area : numpy float32 array
         Pixel area map array

    pathloss : numpy float32 array
         Pathloss correction

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/image.schema"

    def __init__(self, init=None, **kwargs):

        super(ImageModel, self).__init__(init=init, **kwargs)
