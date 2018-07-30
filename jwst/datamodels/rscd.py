from . import model_base


__all__ = ['RSCDModel']


class RSCDModel(model_base.DataModel):
    """
    A data model for the RSCD reference file.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    rscd_table : numpy array
        A table with seven columns, three string-valued that identify which
        row to select, and four float columns containing coefficients.
    """
    schema_url = "rscd.schema.yaml"

    def __init__(self, init=None, rscd_table=None, **kwargs):
        super(RSCDModel, self).__init__(init=init, **kwargs)

        if rscd_table is not None:
            self.rscd_table = rscd_table
