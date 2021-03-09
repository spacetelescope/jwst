from .model_base import JwstDataModel


__all__ = ['ImageModel']


class ImageModel(JwstDataModel):
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

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise

    area : numpy float32 array
         Pixel area map array

    pathloss_point : numpy float32 array
         Pathloss correction for point source

    pathloss_uniform : numpy float32 array
         Pathloss correction for uniform source
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/image.schema"
