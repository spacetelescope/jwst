from __future__ import (
    absolute_import,
    unicode_literals,
    division,
    print_function
)

from . import model_base
from .image import ImageModel


__all__ = ['MultiExposureModel']


class MultiExposureModel(model_base.DataModel):
    """
    A data model for multi-slit images derived from
    numerous exposures. The intent is that all slits
    in this model are of the same source, with each slit
    representing a separate exposure of that source.

    This model has a special member `slits` that can be used to
    deal with an entire slit at a time.  It behaves like a list::

       >>> multislit_model.slits.append(image_model)
       >>> multislit_model.slits[0]
       <ImageModel>

    Also, there is an extra mete attribute, `exposures`. This will
    contain a list of dicts representing the meta attribute from
    the exposures from which each slit has been taken.

    See the module `exp_to_source` for the initial creation of these
    models. This is part of the Level 3 processing of multi-objection
    observations.
    """
    schema_url = "multiexposure.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, ImageModel):
            super(MultiExposureModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.exposures.append(self.exposures.item())
            self.exposures[0].data = init.data
            self.exposures[0].dq = init.dq
            self.exposures[0].err = init.err
            self.exposures[0].relsens = init.relsens
            self.exposures[0].area = init.area
            return

        super(MultiExposureModel, self).__init__(init=init, **kwargs)
