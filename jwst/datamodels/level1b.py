import numpy as np

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

    int_times : table
        The int_times table

    """
    schema_url = "level1b.schema.yaml"

    def __init__(self, init=None, data=None, refout=None, zeroframe=None,
                 group=None, int_times=None, **kwargs):
        super(Level1bModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if refout is not None:
            self.refout = refout

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if group is not None:
            self.group = group

        if int_times is not None:
            self.int_times = int_times

        # zeroframe is a lower dimensional array than
        # the science data. However, its dimensions are not
        # consecutive with data, so the default model
        # creates a wrongly shaped array. If data is given
        # use the appropriate dimensions.
        #
        # TODO: Hacky. Need solution which involves schema
        # specification and embedded in DataModel.
        #if 'zeroframe' not in self.instance and \
        #   'data' in self.instance and \
        #   len(self.data.shape) == 4:
        #    nints, ngroups, ny, nx = self.data.shape
        #    self.zeroframe = np.zeros((nints, ny, nx))
