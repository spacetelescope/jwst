from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['FGSModel']


class FGSModel(model_base.DataModel):
    """
    A data model for input for the FGS pipeline.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data. 4-D

    dq : numpy array
        The data quality array. 2-D.

    plan_star : table
        The planned reference star table

    flight_star : table
        The flight reference star table

    pointing : table
        The pointing table

    centroid : table
        The centroid packet table

    track_sub : table
        The track subarray table

    """
    schema_url = "fgs.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, 
                 plan_star_table=None, flight_star_table=None, 
                 pointing_table=None, centroid_table=None, 
                 track_sub_table=None, **kwargs):
        super(FGSModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if plan_star_table is not None:
            self.plan_star_table = plan_star_table

        if flight_star_table is not None:
            self.flight_star_table = flight_star_table

        if pointing_table is not None:
            self.pointing_table = pointing_table

        if centroid_table is not None:
            self.centroid_table = centroid_table

        if track_sub_table is not None:
            self.track_sub_table = track_sub_table

        # Implicitly create arrays
        self.dq = self.dq

