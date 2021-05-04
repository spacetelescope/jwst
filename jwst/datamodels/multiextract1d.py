from .reference import ReferenceFileModel
from .extract1dimage import Extract1dImageModel


__all__ = ['MultiExtract1dImageModel']


class MultiExtract1dImageModel(ReferenceFileModel):
    """
    A data model for extract_1d reference images.

    This model has a special member `images` that can be used to
    deal with each image separately.  It behaves like a list::

       >>> from jwst.datamodels import Extract1dImageModel
       >>> multiextr1d_img_model = MultiExtract1dImageModel()
       >>> multiextr1d_img_model.images.append(Extract1dImageModel())
       >>> multiextr1d_img_model.images[0]  # doctest: +SKIP
       <Extract1dImageModelModel>

    If `init` is a file name or an `Extract1dImageModel` instance, an empty
    `Extract1dImageModel` will be created and assigned to attribute
    `images[0]`, and the `data` attribute from the input array or
    `Extract1dImageModel` will be copied to the first element of `images`.

    Parameters
    __________
    images.items.data : numpy float32 array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/multiextract1d.schema"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, Extract1dImageModel):
            super(MultiExtract1dImageModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.images.append(self.images.item())
            self.images[0].data = init.data
            return

        super(MultiExtract1dImageModel, self).__init__(init=init, **kwargs)
