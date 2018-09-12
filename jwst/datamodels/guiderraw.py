from .image import ImageModel

__all__ = ['GuiderRawModel']


class GuiderRawModel(ImageModel):
    """
    A data model for FGS pipeline input files

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data. 4-D

    dq: numpy array
        The data quality array. 2-D.

    err: numpy array
        The error array. 4-D.

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

    schema_url = "guider_raw.schema.yaml"
