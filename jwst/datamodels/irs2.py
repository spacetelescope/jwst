from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['IRS2Model']


class IRS2Model(model_base.DataModel):
    """
    A data model for the IRS2 reference file.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    irs2_table : numpy array
        A table with 8 columns and 1458176 (2048 * 712) rows; all values
        are complex.  There are four columns for ALPHA and four for BETA.
    """
    schema_url = "irs2.schema.yaml"

    def __init__(self, init=None, irs2_table=None, **kwargs):
        super(IRS2Model, self).__init__(init=init, **kwargs)

        if irs2_table is not None:
            self.irs2_table = irs2_table
