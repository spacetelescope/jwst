from .model_base import DataModel


__all__ = ['RampModel']


class RampModel(DataModel):
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

    def __init__(self, init=None, **kwargs):
        super(RampModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.pixeldq = self.pixeldq
        self.groupdq = self.groupdq
        self.err = self.err
