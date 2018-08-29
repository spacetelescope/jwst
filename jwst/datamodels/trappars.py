from .reference import ReferenceFileModel


__all__ = ['TrapParsModel']


class TrapParsModel(ReferenceFileModel):
    """
    A data model for trap capture and decay parameters.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    trappars_table : numpy array
        A table with three columns for trap-capture parameters and one
        column for the trap-decay parameter.  Each row of the table is
        for a different trap family.
    """
    schema_url = "trappars.schema.yaml"

    def __init__(self, init=None, trappars_table=None, **kwargs):
        super(TrapParsModel, self).__init__(init=init, **kwargs)

        if trappars_table is not None:
            self.trappars_table = trappars_table
