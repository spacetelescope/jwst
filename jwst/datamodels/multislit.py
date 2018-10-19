from . import model_base
from .image import ImageModel
from .slit import SlitModel, SlitDataModel


__all__ = ['MultiSlitModel']



class MultiSlitModel(model_base.DataModel):
    """
    A data model for multi-slit images.

    This model has a special member `slits` that can be used to
    deal with an entire slit at a time.  It behaves like a list::

       >>> from .slit import SlitModel
       >>> multislit_model = MultiSlitModel()
       >>> multislit_model.slits.append(SlitModel())
       >>> multislit_model.slits[0]
       >>> multislit[0]
       <SlitModel>

    If ``init`` is a file name or an ``ImageModel`` or a ``SlitModel``instance,
    an empty ``SlitModel`` will be created and assigned to attribute ``slits[0]``,
    and the `data`, ``dq``, ``err``, ``var_rnoise``, ``var_poisson``and
    ``relsens`` attributes from the input file or model will be copied to the
    first element of ``slits``.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.
    """
    schema_url = "multislit.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, (SlitModel, ImageModel)):
            super(MultiSlitModel, self).__init__(init=None, **kwargs)
            self.update(init)
            slitdata = SlitDataModel(init)
            self.slits.append(slitdata)
            return

        super(MultiSlitModel, self).__init__(init=init, **kwargs)

    def __getitem__(self, key):
        """
        Returns a metadata value using a dotted name or
        a ``SlitModel``.
        """
        if isinstance(key, str) and key.split('.')[0] == 'meta':
            res = super(MultiSlitModel, self).__getitem__(key)
            return res
        elif isinstance(key, int):
            # Return an instance of a SlitModel
            slit = self.slits[key]  # returns an ObjectNode instance

            kwargs = {}
            items = dict(slit.items())
            for key in items:
                if not key.startswith(('meta', 'extra_fits')):
                    kwargs[key] = items[key]
            s = SlitModel(**kwargs)
            s.update(self)

            if slit.meta.hasattr('wcs'):
                s.meta.wcs = slit.meta.wcs
            return s
        else:
            raise ValueError("Invalid key {0}".format(key))
