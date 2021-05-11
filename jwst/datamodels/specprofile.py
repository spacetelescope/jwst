from .reference import ReferenceFileModel


__all__ = ['SpecProfileModel', 'SpecProfileSingleModel']


class SpecProfileModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS spectral profile reference files.

    This model has a special member `profile` that can be used to
    deal with an entire spectral profile at a time.  It behaves like a list::

       >>> from jwst.datamodels import SpecProfileSingleModel
       >>> specprofile_model = SpecProfileModel()
       >>> specprofile_model.profile.append(SpecProfileSingleModel())
       >>> specprofile_model.profile[0] # doctest: +SKIP
       <SpecProfileSingleModel>

    If `init` is a `SpecProfileSingleModel` instance, an empty `SpecProfileSingleModel`
    will be created and assigned to attribute `profile[0]`, and the `data`
    attribute from the input `SpecProfileSingleModel` instance will be copied to
    the first element of `profile`.  `SpecProfileSingleModel` objects can be appended
    to the `profile` attribute by using its `append` method.

    Parameters
    __________
    profile.items.data : numpy array
         Spectral profile data

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/specprofile.schema"

    def __init__(self, init=None, **kwargs):

        if isinstance(init, SpecProfileSingleModel):
            super(SpecProfileModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.profile.append(self.profile.item())
            self.profile[0].data = init.data

        super(SpecProfileModel, self).__init__(init=init, **kwargs)


class SpecProfileSingleModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS spectral profile data.

    Parameters
    __________
    data : numpy float32 array
         Spectral profile values
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/specprofilesingle.schema"
