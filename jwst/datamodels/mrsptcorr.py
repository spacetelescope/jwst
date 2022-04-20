from .reference import ReferenceFileModel

__all__ = ['MirMrsPtCorrModel']


class MirMrsPtCorrModel(ReferenceFileModel):
    """
    A data model for MIRI mrs IFU across-slice corrections file.

    Parameters
    __________
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        An array-like object containing the pixel-by-pixel spectral leak values
        in units of (MJy / sr) / (DN / sec).

    err : numpy array
        An array-like object containing the uncertainties in the spectral leak
        values, in the same units as the data array.

    dq : numpy array
        An array-like object containing bit-encoded data quality flags,
        indicating problem conditions for values in the data array.

    dq_def : numpy array
        A table-like object containing the data quality definitions table.

    tracor_table : numpy table
         IFU across slice transmission correction

    wavcorr_optical_table : numpy table
         IFU across slice wavelength offset table 1

    wavcorr_xslice_table : numpy table
         IFU across slice wavelength offset table 2

    wavcorr_shift_table : numpy table
         IFU across slice wavelength offset table 3

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/miri_mrsptcorr.schema"
