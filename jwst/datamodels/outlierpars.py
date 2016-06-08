from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['OutlierParsModel']


class OutlierParsModel(model_base.DataModel):
    """
    A data model for outlier detection parameters reference tables.
    """
    schema_url = "outlierpars.schema.yaml"

    def __init__(self, init=None, outlierpars_table=None, **kwargs):
        super(OutlierParsModel, self).__init__(init=init, **kwargs)

        if outlierpars_table is not None:
            self.outlierpars_table = outlierpars_table


class NircamOutlierParsModel(OutlierParsModel):
    """
    A data model for NIRCam outlierpars reference files.
    """
    schema_url = "nircam_outlierpars.schema.yaml"

    def __init__(self, init=None, outlierpars_table=None, **kwargs):
        super(NircamOutlierParsModel, self).__init__(init=init, **kwargs)

        if outlierpars_table is not None:
            self.outlierpars_table = outlierpars_table


class MiriImgOutlierParsModel(OutlierParsModel):
    """
    A data model for MIRI imaging outlierpars reference files.
    """
    schema_url = "mirimg_outlierpars.schema.yaml"

    def __init__(self, init=None, outlierpars_table=None, **kwargs):
        super(MiriImgOutlierParsModel, self).__init__(init=init, **kwargs)

        if outlierpars_table is not None:
            self.outlierpars_table = outlierpars_table
