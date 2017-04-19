from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['TrapsModel']


class TrapsModel(model_base.DataModel):
    """
    A data model for trap capture and decay.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    traps_table : numpy array
        A table with three columns for trap-capture parameters and one
        column for the trap-decay timescale.  Each row of the table is
        for a different trap family.
    """
    schema_url = "traps.schema.yaml"

    def __init__(self, init=None, traps_table=None, **kwargs):
        super(TrapsModel, self).__init__(init=init, **kwargs)

        if traps_table is not None:
            self.traps_table = traps_table
