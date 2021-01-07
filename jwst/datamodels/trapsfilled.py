from .model_base import JwstDataModel


__all__ = ['TrapsFilledModel']


class TrapsFilledModel(JwstDataModel):
    """
    A data model for the number of traps filled for a detector, for
    persistence.

    Parameters
    __________
    data : numpy float32 array
        Traps filled
        The map of the number of traps filled over the detector, with
        one plane for each "trap family."
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/trapsfilled.schema"
