from .model_base import DataModel


__all__ = ['ImageModel']


class ImageModel(DataModel):
    """
    A data model for 2D images.

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

    pathloss_pointsource : numpy float32 array
         Pathloss correction

    relsens : numpy table
         relative sensitivity table

    relsens2d : numpy float32 array
         Sensitivity array

    var_poisson : numpy float32 array
         variance due to poisson noise

    var_rnoise : numpy float32 array
         variance due to read noise
    """
    schema_url = "image.schema.yaml"

    def __init__(self, init=None, **kwargs):

        super(ImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
