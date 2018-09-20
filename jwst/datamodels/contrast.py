from .model_base import DataModel

__all__ = ['ContrastModel']


class ContrastModel(DataModel):
    """
    A data model for coronagraphic contrast curve files.
    """
    schema_url = "contrast.schema.yaml"
