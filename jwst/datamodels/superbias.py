from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask


__all__ = ['SuperBiasModel']


class SuperBiasModel(ReferenceFileModel):
    """
    A data model for 2D super-bias images.
    """
    schema_url = "superbias.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(SuperBiasModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
