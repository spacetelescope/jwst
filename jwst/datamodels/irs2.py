from .model_base import DataModel


__all__ = ['IRS2Model']


class IRS2Model(DataModel):
    """
    A data model for the IRS2 refpix reference file.

    Parameters
    __________
    irs2_table : numpy table
         Reference file for IRS2 refpix correction
         A table with 8 columns and 2916352 (2048 * 712 * 2) rows.  All
         values are float, but these are interpreted as alternating real
         and imaginary parts (real, imag, real, imag, ...) of complex
         values.  There are four columns for ALPHA and four for BETA.
    """
    schema_url = "irs2.schema"
