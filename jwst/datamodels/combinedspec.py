from .model_base import DataModel

__all__ = ['CombinedSpecModel']


class CombinedSpecModel(DataModel):
    """
    A data model for combined 1D spectra.
    """
    schema_url = "combinedspec.schema.yaml"
