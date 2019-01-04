from .model_base import DataModel

__all__ = ['Level1bModel']


class Level1bModel(DataModel):
    """
    A data model for raw 4D ramps level-1b products.

    Parameters
    __________
    data : numpy uint16 array
         The science data

    zeroframe : numpy uint16 array
         Zeroframe array

    refout : numpy uint16 array
         Reference Output

    group : numpy table
         group parameters table

    int_times : numpy table
         table of times for each integration

    """
    schema_url = "level1b.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(Level1bModel, self).__init__(init=init, **kwargs)

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
