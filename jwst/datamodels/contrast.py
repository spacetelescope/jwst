from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base

__all__ = ['ContrastModel']


class ContrastModel(model_base.DataModel):
    """
    A data model for coronagraphic contrast curve files.
    """
    schema_url = "contrast.schema.yaml"

    def __init__(self, init=None, contrast_table=None, **kwargs):
        super(ContrastModel, self).__init__(init=init, **kwargs)

        if contrast_table is not None:
            self.contrast_table = contrast_table
