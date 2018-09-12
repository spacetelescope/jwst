from .model_base import DataModel


__all__ = ['ImageModel']


class ImageModel(DataModel):
    """
    A data model for 2D images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    kwargs: numpy array
        The name and value of any array mentioned in the schema to be
        initialized through the function call.
    """
    schema_url = "image.schema.yaml"

    def __init__(self, init=None, **kwargs):

        super(ImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

        if self.hasattr('dq') and self.hasattr('dq_def'):
            self.dq = dynamic_mask(self)

