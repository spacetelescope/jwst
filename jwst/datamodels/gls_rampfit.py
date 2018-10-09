from .model_base import DataModel

__all__ = ['GLS_RampFitModel']


class GLS_RampFitModel(DataModel):
    """
    A data model for the optional output of the ramp fitting step
    for the GLS algorithm.
    """
    schema_url = "gls_rampfit.schema.yaml"
