from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base
from .spec import SpecModel


__all__ = ['MultiSpecModel']


class MultiSpecModel(model_base.DataModel):
    """
    A data model for multi-spec images.

    This model has a special member `spec` that can be used to
    deal with an entire spectrum at a time.  It behaves like a list::

       >>> multispec_model.spec.append(spec_model)
       >>> multispec_model.spec[0]
       <SpecModel>

    If `init` is a `SpecModel` instance, an empty `SpecModel` will be
    created and assigned to attribute `spec[0]`, and the `spec_table`
    attribute from the input `SpecModel` instance will be copied to
    the first element of `spec`.  `SpecModel` objects can be appended
    to the `spec` attribute by using its `append` method.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    Examples
    --------
    >>> output_model = models.MultiSpecModel()
    >>> spec = models.SpecModel()           # for the default data type
    >>> for slit in input_model.slits:
    >>>     slitname = slit.name
    >>>     slitmodel = ExtractModel()
    >>>     slitmodel.fromJSONFile(extref, slitname)
    >>>     column, wavelength, countrate = slitmodel.extract(slit.data)
    >>>     otab = np.array(zip(column, wavelength, countrate),
    >>>                     dtype=spec.spec_table.dtype)
    >>>     spec = models.SpecModel(spec_table=otab)
    >>>     output_model.spec.append(spec)
    """
    schema_url = "multispec.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, SpecModel):
            super(MultiSpecModel, self).__init__(init=None, **kwargs)
            self.spec.append(self.spec.item())
            self.spec[0].spec_table = init.spec_table
            return

        super(MultiSpecModel, self).__init__(init=init, **kwargs)
