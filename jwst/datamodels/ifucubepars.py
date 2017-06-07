from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['IFUCubeParsModel']


class IFUCubeParsModel(model_base.DataModel):
    """
    A data model for IFU Cube  parameters reference tables.
    """
    schema_url = "ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None, **kwargs):
        super(IFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table


class NirspecIFUCubeParsModel(IFUCubeParsModel):
    """
    A data model for Nirspec ifucubepars reference files.
    """
    schema_url = "nirspec_ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None, **kwargs):
        super(NirspecIFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table


class MiriIFUCubeParsModel(IFUCubeParsModel):
    """
    A data model for MIRI mrs ifucubepars reference files.
    """
    schema_url = "miri_ifucubepars.schema.yaml"

    def __init__(self, init=None, ifucubepars_table=None, **kwargs):
        super(MiriIFUCubeParsModel, self).__init__(init=init, **kwargs)

        if ifucubepars_table is not None:
            self.ifucubepars_table = ifucubepars_table
