from .model_base import JwstDataModel
from .image import ImageModel
from .slit import SlitModel, SlitDataModel


__all__ = ['MultiSlitModel']


class MultiSlitModel(JwstDataModel):
    """
    A data model for multi-slit images.

    This model has a special member `slits` that can be used to
    deal with an entire slit at a time.  It behaves like a list::

       >>> from jwst.datamodels import SlitModel
       >>> multislit_model = MultiSlitModel()
       >>> multislit_model.slits.append(SlitModel())
       >>> multislit_model[0]
       <SlitModel>

    If ``init`` is a file name or an ``ImageModel`` or a ``SlitModel``instance,
    an empty ``SlitModel`` will be created and assigned to attribute ``slits[0]``,
    and the `data`, ``dq``, ``err``, ``var_rnoise``, and ``var_poisson``
    attributes from the input file or model will be copied to the
    first element of ``slits``.

    Parameters
    __________
    slits.items.data : numpy float32 array
         The science data

    slits.items.dq : numpy uint32 array
         Data quality array

    slits.items.err : numpy float32 array
         Error array

    slits.items.var_poisson : numpy float32 array
         variance due to poisson noise

    slits.items.var_rnoise : numpy float32 array
         variance due to read noise

    slits.items.wavelength : numpy float32 array
         Wavelength array, corrected for zero-point

    slits.items.barshadow : numpy float32 array
         Bar shadow correction

    slits.items.flatfield_point : numpy float32 array
         flatfield array for point source

    slits.items.flatfield_uniform : numpy float32 array
         flatfield array for uniform source

    slits.items.pathloss_point : numpy float32 array
         pathloss array for point source

    slits.items.pathloss_uniform : numpy float32 array
         pathloss array for uniform source

    slits.items.photom_point : numpy float32 array
         photom array for point source

    slits.items.photom_uniform : numpy float32 array
         photom array for uniform source

    slits.items.area : numpy float32 array
         Pixel area map array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/multislit.schema"

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
            # This only executes when the top object level
            # is called: object[1].key not object.slits[key]
            try:
                slit = self.slits[key]  # returns an ObjectNode instance
            except IndexError:
                raise ("Slit {0} doesn't exist".format(key))
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
