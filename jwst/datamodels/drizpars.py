from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['DrizParsModel']


class DrizParsModel(model_base.DataModel):
    """
    A data model for drizzle parameters reference tables.
    """
    schema_url = "drizpars.schema.yaml"

    def __init__(self, init=None, drizpars_table=None, **kwargs):
        super(DrizParsModel, self).__init__(init=init, **kwargs)

        if drizpars_table is not None:
            self.drizpars_table = drizpars_table


class NircamDrizParsModel(DrizParsModel):
    """
    A data model for NIRCam drizpars reference files.
    """
    schema_url = "nircam_drizpars.schema.yaml"

    def __init__(self, init=None, drizpars_table=None, **kwargs):
        super(NircamDrizParsModel, self).__init__(init=init, **kwargs)

        if drizpars_table is not None:
            self.drizpars_table = drizpars_table


class MiriImgDrizParsModel(DrizParsModel):
    """
    A data model for MIRI imaging drizpars reference files.
    """
    schema_url = "mirimg_drizpars.schema.yaml"

    def __init__(self, init=None, drizpars_table=None, **kwargs):
        super(MiriImgDrizParsModel, self).__init__(init=init, **kwargs)

        if drizpars_table is not None:
            self.drizpars_table = drizpars_table
