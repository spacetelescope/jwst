from .model_base import JwstDataModel


__all__ = ['IRS2Model']


class IRS2Model(JwstDataModel):
    """
    A data model for the IRS2 refpix reference file.

    Parameters
    __________
    irs2_table : numpy table
        Table for IRS2 refpix correction.
        A table with 8 columns and 2916352 (2048 * 712 * 2) rows.  All
        values are float, but these are interpreted as alternating real
        and imaginary parts (real, imag, real, imag, ...) of complex
        values.  There are four columns for ALPHA and four for BETA.

    dq_table : data quality info table
        Table for identifying bad reference pixels.
        A table with three columns (OUTPUT, ODD_EVEN, and MASK) and
        eight rows.
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/irs2.schema"
