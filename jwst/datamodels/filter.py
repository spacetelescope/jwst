from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['FilterModel']


class FilterModel(model_base.DataModel):
    """
    A data model for filter throughput.
    """
    schema_url = "filter.schema.yaml"

    def __init__(self, init=None, filter_table=None, **kwargs):
        super(FilterModel, self).__init__(init=init, **kwargs)

        if filter_table is not None:
            self.filter_table = filter_table
