from . import model_base
from .dynamicdq import dynamic_mask


__all__ = ['PathlossModel']


class PathlossModel(model_base.DataModel):
    """
    A data model for pathloss correction information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    pointsource : numpy array
        Array defining the pathloss parameter for point sources.

    psvar : numpy array
        Variance array.

    uniform : numpy array
        Pathloss parameter for uniform illumination
    """
    schema_url = "pathloss.schema.yaml"

    def __init__(self, init=None, pointsource=None, psvar=None, uniform=None,
                 **kwargs):
        super(PathlossModel, self).__init__(init=init, **kwargs)

        if pointsource is not None:
            self.pointsource = pointsource

        if psvar is not None:
            self.psvar = psvar

        if uniform is not None:
            self.uniform = uniform

        # Implicitly create arrays
