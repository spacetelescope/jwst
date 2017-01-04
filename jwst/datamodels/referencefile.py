from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['ReferencefileModel']


class ReferencefileModel(model_base.DataModel):
    """
    A data model for 2D images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.
    """
    schema_url = "referencefile.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(ReferencefileModel, self).__init__(init=init, **kwargs)

