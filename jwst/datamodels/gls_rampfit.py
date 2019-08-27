from .model_base import DataModel

__all__ = ['GLS_RampFitModel']


class GLS_RampFitModel(DataModel):
    """
    A data model for the optional output of the ramp fitting step
    for the GLS algorithm.

    Parameters
    __________
    yint : numpy float32 array
         Y-intercept

    sigyint : numpy float32 array
         Sigma for Y-intercept

    pedestal : numpy float32 array
         Pedestal

    crmag : numpy float32 array
         CR magnitudes

    sigcrmag : numpy float32 array
         Sigma for CR magnitudes
    """
    schema_url = "gls_rampfit.schema"
