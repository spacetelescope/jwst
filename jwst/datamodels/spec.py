from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['SpecModel']


class SpecModel(model_base.DataModel):
    """
    A data model for 1D spectra.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    spec_table : numpy array
        A table with at least four columns:  wavelength, flux, an error
        estimate for the flux, and data quality flags.
    """
    schema_url = "spec.schema.yaml"

    def __init__(self, init=None, spec_table=None, **kwargs):
        super(SpecModel, self).__init__(init=init, **kwargs)

        if spec_table is not None:
            self.spec_table = spec_table
