from .ramp import RampModel


__all__ = ['MIRIRampModel']


class MIRIRampModel(RampModel):
    """
    A data model for MIRI ramps. Includes the ``refout`` array.

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

    refout : numpy array
        The array of reference output data.

    group : table
        The group parameters table.

    """
    schema_url = "miri_ramp.schema.yaml"

