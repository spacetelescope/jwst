from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['PixelAreaModel']


class PixelAreaModel(model_base.DataModel):
    """
    A data model for the pixel area map
    """
    schema_url = "pixelarea.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(PixelAreaModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data
