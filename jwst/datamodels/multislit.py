from . import model_base
from .image import ImageModel
from .slit import SlitModel, SlitDataModel


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
        Any of the initializers supported by `~jwst.datamodels.DataModel`.
    """
    schema_url = "multislit.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, SlitModel):
            super(MultiSlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.slits.append(self.slits.item())
            self.slits[0].data = init.data
            self.slits[0].dq = init.dq
            self.slits[0].err = init.err
            self.slits[0].relsens = init.relsens
            self.slits[0].area = init.area
            self.slits[0].wavelength = init.wavelength
            return

        super(MultiSlitModel, self).__init__(init=init, **kwargs)

    def __getitem__(self, key):
        """
        Get a metadata value using a dotted name.
        """
        def _get_data_keys(*obj):
            return obj

        if isinstance(key, six.string_types) and key.split('.') == 'meta':
            super(MultiSlitModel, self).__getitem__(key)
        elif isinstance(key, int):
            # Return an instance of a SlitModel
            data_keys = _get_data_keys(*self.slits[key])
            kwargs = dict(((k, getattr(self.slits[key], k)) for k in data_keys))
            s = SlitModel(**kwargs)
            s.update(self)#, only='PRIMARY')
            return s
        else:
            raise ValueError("Invalid key {0}".format(key))
