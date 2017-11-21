from . import model_base

__all__ = ['BarshadowModel']


class BarshadowModel(model_base.DataModel):
    """
    A data model for Bar Shadow correction information.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        Array defining the bar shadow correction as a function of Y and
        wavelength.

    variance : numpy array
        Variance array.

    """
    schema_url = "barshadow.schema.yaml"

    def __init__(self, init=None, data1x1=None, var1x1=None, data1x3=None,
                 var1x3=None, **kwargs):
        super(BarshadowModel, self).__init__(init=init, **kwargs)

        if data1x1 is not None:
            self.data1x1 = data1x1

        if var1x1 is not None:
            self.var1x1 = var1x1

        if data1x3 is not None:
            self.data1x3 = data1x3

        if var1x3 is not None:
            self.var1x3 = var1x3

# Implicitly create arrays
