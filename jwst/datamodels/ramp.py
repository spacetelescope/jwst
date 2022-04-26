import numpy as np

from .model_base import JwstDataModel
from ..lib.basic_utils import deprecate_class

__all__ = ['RampModel']


class RampModel(JwstDataModel):
    """
    A data model for 4D ramps.

    Parameters
    __________
    data : numpy float32 array
         The science data

    pixeldq : numpy uint32 array
         2-D data quality array for all planes

    groupdq : numpy uint8 array
         4-D data quality array for each plane

    err : numpy float32 array
         Error array

    zeroframe : numpy float32 array
         Zeroframe array

    group : numpy table
         group parameters table

    int_times : numpy table
         table of times for each integration

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/ramp.schema"

    def __init__(self, init=None, **kwargs):
        super(RampModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.pixeldq = self.pixeldq
        self.groupdq = self.groupdq
        self.err = self.err

        if isinstance(init, tuple) or self.meta.exposure.zero_frame is True:
            try:
                self.getarray_noinit("zeroframe")
            except AttributeError:
                # If "zeroframe" is not in the instance, create a zero array with
                # the correct dimensions.
                nints, ngroups, nrows, ncols = self.data.shape
                dims = (nints, nrows, ncols)
                self.zeroframe = np.zeros(dims, dtype=self.data.dtype)


@deprecate_class(RampModel)
class MIRIRampModel:
    """A data model for 4D MIRI ramps.

    This model has been deprecated. Please use `RampModel` instead.
    """
