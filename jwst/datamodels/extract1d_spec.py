from .reference import ReferenceFileModel

__all__ = ['Extract1dIFUModel']


class Extract1dIFUModel(ReferenceFileModel):
    """
    A data model for IFU MIRI and NIRSpec extract 1d reference files.

    Parameters
    __________
    extract1d_params : numpy table
        Basic extract 1D parameters
        - region_type: ascii
        - subtract_background: bool
        - method: ascii
        - subpixels: int16
    extract1d_table : numpy table
        extract1d parameters for varying wavelengths
        A table-like object containing extract 1d parameters
        based on wavelength
        - wavelength: float32 1D array
        - radius: float32 1D array
        - inner_bkg: float32 1D array
        - outer_bkg: float32 1D array
        - axis_ratio: float32 1D array
        - axis_pa: float32 1D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/extract1difu.schema"
