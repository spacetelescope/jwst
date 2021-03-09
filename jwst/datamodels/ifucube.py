from .model_base import JwstDataModel


__all__ = ['IFUCubeModel']


class IFUCubeModel(JwstDataModel):
    """
    A data model for 3D IFU  cubes.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array

    weightmap : numpy float32 array
         Weight map of coverage

    wavetable : numpy table
         Wavelength value for slices
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ifucube.schema"

    def __init__(self, init=None, **kwargs):
        super(IFUCubeModel, self).__init__(init=init, **kwargs)

       # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
