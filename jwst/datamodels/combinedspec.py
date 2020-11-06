from .model_base import JwstDataModel


__all__ = ['CombinedSpecModel']


class CombinedSpecModel(JwstDataModel):
    """
    A data model for combined 1D spectra.

    Parameters
    __________
    spec_table : numpy table
         Combined, extracted spectral data table
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/combinedspec.schema"
