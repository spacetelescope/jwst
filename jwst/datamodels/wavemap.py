from .reference import ReferenceFileModel


__all__ = ['WaveMapModel', 'WaveMapSingleModel']


class WaveMapModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS wavelength map reference files.

    This model has a special member `map` that can be used to
    deal with an entire wavelength map at a time.  It behaves like a list::

       >>> from jwst.datamodels import WaveMapSingleModel, WaveMapModel
       >>> wavemap_model = WaveMapModel()
       >>> wavemap_model.map.append(WaveMapSingleModel())
       >>> wavemap_model.map[0] # doctest: +SKIP
       <WaveMapSingleModel>

    If `init` is a `WaveMapSingleModel` instance, an empty `WaveMapSingleModel`
    will be created and assigned to attribute `map[0]`, and the `data`
    attribute from the input `WaveMapSingleModel` instance will be copied to
    the first element of `map`.  `WaveMapSingleModel` objects can be appended
    to the `map` attribute by using its `append` method.

    Parameters
    __________
    map.items.data : numpy data array
         Wavelength map data

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/wavemap.schema"

    def __init__(self, init=None, **kwargs):

        if isinstance(init, WaveMapSingleModel):
            super(WaveMapModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.map.append(self.map.item())
            self.map[0].data = init.data

        super(WaveMapModel, self).__init__(init=init, **kwargs)


class WaveMapSingleModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS wavelength map data.

    Parameters
    __________
    data : numpy float32 array
         Wavelength values
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/wavemapsingle.schema"
