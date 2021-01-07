from .model_base import JwstDataModel
from .spec import SpecModel


__all__ = ['MultiSpecModel']


class MultiSpecModel(JwstDataModel):
    """
    A data model for multi-spec images.

    This model has a special member `spec` that can be used to
    deal with an entire spectrum at a time.  It behaves like a list::

       >>> from . import SpecModel
       >>> multispec_model = MultiSpecModel()
       >>> multispec_model.spec.append(SpecModel())
       >>> multispec_model.spec[0] # doctest: +SKIP
       <SpecModel>

    If `init` is a `SpecModel` instance, an empty `SpecModel` will be
    created and assigned to attribute `spec[0]`, and the `spec_table`
    attribute from the input `SpecModel` instance will be copied to
    the first element of `spec`.  `SpecModel` objects can be appended
    to the `spec` attribute by using its `append` method.

    Parameters
    __________
    int_times : numpy table
         table of times for each integration

    spec.items.spec_table : numpy table
         Extracted spectral data table

    Examples
    --------
    >>> output_model = MultiSpecModel()
    >>> spec = SpecModel()       # for the default data type
    >>> for slit in input_model.slits:  # doctest: +SKIP
    ...     slitname = slit.name
    ...     slitmodel = ExtractModel()
    ...     slitmodel.fromJSONFile(extref, slitname)
    ...     column, wavelength, countrate = slitmodel.extract(slit.data)
    ...     otab = np.array(zip(column, wavelength, countrate),
    ...                     dtype=spec.spec_table.dtype)
    ...     spec = datamodels.SpecModel(spec_table=otab)
    ...     output_model.spec.append(spec)
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/multispec.schema"

    def __init__(self, init=None, **kwargs):

        if isinstance(init, SpecModel):
            super(MultiSpecModel, self).__init__(init=None, **kwargs)
            self.spec.append(self.spec.item())
            self.spec[0].spec_table = init.spec_table
            return

        super(MultiSpecModel, self).__init__(init=init, **kwargs)
