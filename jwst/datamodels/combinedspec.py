from .model_base import DataModel

__all__ = ['CombinedSpecModel']


class CombinedSpecModel(DataModel):
    """
    A data model for combined 1D spectra.

    Attributes
    __________
    spec_table : numpy table
         Combined, extracted spectral data table
    """
    schema_url = "combinedspec.schema.yaml"
