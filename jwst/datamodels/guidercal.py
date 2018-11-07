from .model_base import DataModel

__all__ = ['GuiderCalModel']


class GuiderCalModel(DataModel):
    """
    A data model for FGS pipeline output files

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

    schema_url = "guider_cal.schema.yaml"

    def __init__(self, init=None, **kwargs):

        super(GuiderCalModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
