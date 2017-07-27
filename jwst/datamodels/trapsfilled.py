from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['TrapsFilledModel']

class TrapsFilledModel(model_base.DataModel):
    """
    A data model for the number of traps filled for a detector, for
    persistence.

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The map of the number of traps filled over the detector, with
        one plane for each "trap family."
    """
    schema_url = "trapsfilled.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(TrapsFilledModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
