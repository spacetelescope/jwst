from . import model_base


__all__ = ['RampModel']


class RampModel(model_base.DataModel):
    """
    A data model for 4D ramps.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    pixeldq : numpy array
        2-D data quality array.

    groupdq : numpy array
        3-D or 4-D data quality array.

    err : numpy array
        The error array.

    group : table
        The group parameters table

    int_times : table
        The int_times table

    """
    schema_url = "ramp.schema.yaml"

    def __init__(self, init=None, data=None, pixeldq=None, groupdq=None,
                 err=None, zeroframe=None, group=None, int_times=None,
                 **kwargs):
        super(RampModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if pixeldq is not None:
            self.pixeldq = pixeldq

        if groupdq is not None:
            self.groupdq = groupdq

        if err is not None:
            self.err = err

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if group is not None:
            self.group = group

        if int_times is not None:
            self.int_times = int_times

        # Implicitly create arrays
        self.pixeldq = self.pixeldq
        self.groupdq = self.groupdq
        self.err = self.err
