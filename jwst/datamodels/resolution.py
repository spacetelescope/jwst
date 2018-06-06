from .reference import ReferenceFileModel

__all__ = ['ResolutionModel', 'MiriResolutionModel']


class ResolutionModel(ReferenceFileModel):
    """
    A data model for Spectral Resolution  parameters reference tables.
    """
    schema_url = "resolution.schema.yaml"

    def __init__(self, init=None, resolution_table=None, **kwargs):
        super(ResolutionModel, self).__init__(init=init, **kwargs)

        if resolution_table is not None:
            self.resolution_table = ifucubepars_table


class MiriResolutionModel(ResolutionModel):
    """
    A data model for MIRI Resolution reference files.

    Parameters
    ----------
    init : any
       Any of the initializers supported by '~jwst.datamodels.DataModel'

    resolving_power_table : table
      A table containing resolving power of the MRS. THe table consist of 11
      columns and 12 rows. Each row corresponds to a band. The columns give the
      name of band, central wavelength, and polynomial coefficeints (a,b,c)
      needed to obtain the limits and average value of the spectral resolution.

    psf_fwhm_alpha_table : table
      A table with 5 columns. Column 1 gives the cutoff wavelength where the
      polynomials describing alpha FWHM change. Columns 2 and 3 give the
      polynomial cofficients (a,b) describing alpha FWHM for wavelengths
      shorter than cuttoff. Columns 4 and 5 give the polynomial
      coefficients (a,b) describing alpha FWHM for wavelengths longer than the
      cutoff.

    psf_fwhm_beta_table : table
      A table with 5 columns. Column 1 gives the cutoff wavelength where the
      polynomials describing alpha FWHM change. Columns 2 and 3 give the
      polynomial cofficients (a,b) describing beta FWHM for wavelengths shorter
      than cuttoff. Columns 4 and 5 give the polynomial coefficients (a,b)
      describing beta FWHM for wavelengths longer than the cutoff.
    """
    schema_url = "miri_resolution.schema.yaml"

    def __init__(self, init=None, resolving_power_table=None, psf_fwhm_alpha_table=None,
                 psf_fwhm_beta_table=None, **kwargs):
        super(MiriResolutionModel, self).__init__(init=init, **kwargs)

        if resolving_power_table is not None:
            self.resolving_power_table = resolving_power_table

        if psf_fwhm_alpha_table is not None:
            self.psf_fwhm_alpha_table = psf_fwhm_alpha_table

        if psf_fwhm_beta_table is not None:
            self.psf_fwhm_beta_table = psf_fwhm_beta_table


