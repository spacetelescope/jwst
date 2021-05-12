from .reference import ReferenceFileModel


__all__ = ['SpecTraceModel', 'SpecTraceSingleModel']


class SpecTraceModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS spectral trace reference files.

    This model has a special member `trace` that can be used to
    deal with an entire spectral trace at a time.  It behaves like a list::

       >>> from jwst.datamodels import SpecTraceSingleModel
       >>> spectrace_model = SpecTraceModel()
       >>> spectrace_model.trace.append(SpecTraceSingleModel())
       >>> spectrace_model.trace[0] # doctest: +SKIP
       <SpecTraceSingleModel>

    If `init` is a `SpecTraceSingleModel` instance, an empty `SpecTraceSingleModel`
    will be created and assigned to attribute `trace[0]`, and the `data`
    attribute from the input `SpecTraceSingleModel` instance will be copied to
    the first element of `trace`.  `SpecTraceSingleModel` objects can be appended
    to the `trace` attribute by using its `append` method.

    Parameters
    __________
    trace.items.data : numpy table
         Spectral trace data

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/spectrace.schema"

    def __init__(self, init=None, **kwargs):

        if isinstance(init, SpecTraceSingleModel):
            super(SpecTraceModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.trace.append(self.trace.item())
            self.trace[0].data = init.data

        super(SpecTraceModel, self).__init__(init=init, **kwargs)


class SpecTraceSingleModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS spectral trace data.

    Parameters
    __________
    data : numpy table
         Trace values
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/spectracesingle.schema"
