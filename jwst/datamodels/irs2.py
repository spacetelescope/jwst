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
        A table with 8 columns and 2916352 (2048 * 712 * 2) rows.  All
        values are float, but these are interpreted as alternating real
        and imaginary parts (real, imag, real, imag, ...) of complex
        values.  There are four columns for ALPHA and four for BETA.
    """
    schema_url = "irs2.schema.yaml"

    def __init__(self, init=None, irs2_table=None, **kwargs):
        super(IRS2Model, self).__init__(init=init, **kwargs)

        if irs2_table is not None:
            self.irs2_table = irs2_table
