from .model_base import DataModel

__all__ = ['TrapsFilledModel']

class TrapsFilledModel(DataModel):
    """
    A data model for the number of traps filled for a detector, for
    persistence.

    Attributes
    __________
    data : numpy float32 array
        Traps filled
        The map of the number of traps filled over the detector, with
        one plane for each "trap family."
    """
    schema_url = "trapsfilled.schema.yaml"
