from .model_base import DataModel

__all__ = ['ContrastModel']


class ContrastModel(DataModel):
    """
    A data model for coronagraphic contrast curve files.

    Parameters
    __________
    contrast_table : numpy table
         Contrast curve table
    """
    schema_url = "contrast.schema.yaml"
