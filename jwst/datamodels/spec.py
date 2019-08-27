from .model_base import DataModel


__all__ = ['SpecModel']


class SpecModel(DataModel):
    """
    A data model for 1D spectra.

    Parameters
    __________
    spec_table : numpy table
        Extracted spectral data table
        A table with at least four columns:  wavelength, flux, an error
        estimate for the flux, and data quality flags.
    """
    schema_url = "spec.schema"
