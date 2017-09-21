from __future__ import absolute_import, unicode_literals, division, print_function

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

    def __init__(self, init=None, data=None, variance=None,
                 **kwargs):
        super(BarshadowModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if variance is not None:
            self.variance = variance

        # Implicitly create arrays
