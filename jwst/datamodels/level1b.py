from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['Level1bModel']


class Level1bModel(model_base.DataModel):
    """
    A data model for raw 4D ramps level-1b products.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data

    zeroframe : numpy array
        The zero-frame data

    refout : numpy array
        The MIRI reference output data

    group : table
        The group parameters table

    """
    schema_url = "level1b.schema.yaml"

    def __init__(self, init=None, data=None, refout=None, zeroframe=None,
                 group=None, **kwargs):
        super(Level1bModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if refout is not None:
            self.refout = refout

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if group is not None:
            self.group = group

