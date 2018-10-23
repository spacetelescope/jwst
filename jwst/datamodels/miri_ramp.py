from .ramp import RampModel


__all__ = ['MIRIRampModel']


class MIRIRampModel(RampModel):
    """
    A data model for MIRI ramps. Includes the ``refout`` array.

    Attributes
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

    refout : numpy float32 array
         Reference Output

    group : numpy table
         group parameters table

    """
    schema_url = "miri_ramp.schema.yaml"

