from .model_base import JwstDataModel


__all__ = ['GuiderRawModel', 'GuiderCalModel']


class GuiderRawModel(JwstDataModel):
    """
    A data model for Guide Star pipeline raw data files

    Parameters
    __________
    data : numpy float32 array
         The science data

    err : numpy float32 array
         Error array

    dq : numpy uint32 array
         Data quality array

    planned_star_table : numpy table
         Planned reference star table

    flight_star_table : numpy table
         Flight reference star table

    pointing_table : numpy table
         Pointing table

    centroid_table : numpy table
         Centroid packet table

    track_sub_table : numpy table
         Track subarray data table
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/guider_raw.schema"

    def __init__(self, init=None, **kwargs):

        super(GuiderRawModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class GuiderCalModel(JwstDataModel):
    """
    A data model for Guide Star pipeline calibrated files

    Parameters
    __________
    data : numpy float32 array
         The science data

    err : numpy float32 array
         Error array

    dq : numpy uint32 array
         Data quality array

    planned_star_table : numpy table
         Planned reference star table

    flight_star_table : numpy table
         Flight reference star table

    pointing_table : numpy table
         Pointing table

    centroid_table : numpy table
         Centroid packet table

    track_sub_table : numpy table
         Track subarray data table
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/guider_cal.schema"

    def __init__(self, init=None, **kwargs):

        super(GuiderCalModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
