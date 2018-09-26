from .model_base import DataModel


__all__ = ['SpecModel']


class SpecModel(DataModel):
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
