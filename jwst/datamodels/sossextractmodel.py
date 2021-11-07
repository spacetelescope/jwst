from .model_base import JwstDataModel

__all__ = ['SossExtractModel']


class SossExtractModel(JwstDataModel):
    """
    A data model to hold NIRISS SOSS extraction model arrays.
    For each order, stores the model trace per integration and
    aperture pixel weights for each order extraction.

    This model is written to explicitly handle each of the three orders.

    Parameters
    __________
    order1 : numpy float32 array
        3-D array of the 2-D model trace for each integration,
        for spectral order 1

    order2 : numpy float32 array
        3-D array of the 2-D model trace for each integration,
        for spectral order 2

    order3 : numpy float32 array
        3-D array of the 2-D model trace for each integration,
        for spectral order 3

    aperture1 : numpy float32 array
        2-D array storing the pixel weights for box-extracting
        spectral order 1

    aperture2 : numpy float32 array
        2-D array storing the pixel weights for box-extracting
        spectral order 2

    aperture3 : numpy float32 array
        2-D array storing the pixel weights for box-extracting
        spectral order 3
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/sossextractmodel.schema"
