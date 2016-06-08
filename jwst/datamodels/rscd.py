from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['RSCD_Model']


class RSCD_Model(model_base.DataModel):
    """
    A data model for the RSCD reference file.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    rscd_table : numpy array
        A table with six columns, two string-valued that identify which
        row to select, and four float columns containing coefficients.
    """
    schema_url = "rscd.schema.yaml"

    def __init__(self, init=None, rscd_table=None, **kwargs):
        super(RSCD_Model, self).__init__(init=init, **kwargs)

        if rscd_table is not None:
            self.rscd_table = rscd_table
