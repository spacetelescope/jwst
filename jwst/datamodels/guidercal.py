from .image import ImageModel

__all__ = ['GuiderCalModel']


class GuiderCalModel(ImageModel):
    """
    A data model for FGS pipeline output files

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data. 3-D

    dq: numpy array
        The data quality array. 2-D

    err: numpy array
        The error array. 3-D

    plan_star_table: table
        The planned reference star table

    flight_star_table: table
        The flight reference star table

    pointing_table: table
        The pointing table

    centroid_table: table
        The centroid packet table

    track_sub_table: table
        The track subarray table
    """

    schema_url = "guider_cal.schema.yaml"
