from .reference import ReferenceImageModel


__all__ = ['SuperBiasModel']


class SuperBiasModel(ReferenceImageModel):
    """
    A data model for 2D super-bias images.
    """
    schema_url = "superbias.schema.yaml"
