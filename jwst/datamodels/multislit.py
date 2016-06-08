from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base
from .image import ImageModel


__all__ = ['MultiSlitModel']


class MultiSlitModel(model_base.DataModel):
    """
    A data model for multi-slit images.

    This model has a special member `slits` that can be used to
    deal with an entire slit at a time.  It behaves like a list::

       >>> multislit_model.slits.append(image_model)
       >>> multislit_model.slits[0]
       <ImageModel>

    If `init` is a file name or an `ImageModel` instance, an empty
    `ImageModel` will be created and assigned to attribute `slits[0]`,
    and the `data`, `dq`, `err`, and `relsens` attributes from the
    input file or `ImageModel` will be copied to the first element of
    `slits`.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.
    """
    schema_url = "multislit.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, ImageModel):
            super(MultiSlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.slits.append(self.slits.item())
            self.slits[0].data = init.data
            self.slits[0].dq = init.dq
            self.slits[0].err = init.err
            self.slits[0].relsens = init.relsens
            self.slits[0].area = init.area
            return

        super(MultiSlitModel, self).__init__(init=init, **kwargs)
