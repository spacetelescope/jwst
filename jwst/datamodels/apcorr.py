from .reference import ReferenceFileModel

__all__ = ['FgsImgApcorrModel', 'MirImgApcorrModel',
           'NrcImgApcorrModel', 'NisImgApcorrModel',
           'MirLrsApcorrModel', 'MirMrsApcorrModel',
           'NrcWfssApcorrModel', 'NisWfssApcorrModel',
           'NrsMosApcorrModel', 'NrsIfuApcorrModel', 'NrsFsApcorrModel']


class FgsImgApcorrModel(ReferenceFileModel):
    """
    A data model for FGS imaging apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - eefraction: float32
        - radius: float32
        - apcorr: float32
        - skyin: float32
        - skyout: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/fgsimg_apcorr.schema"


class MirImgApcorrModel(ReferenceFileModel):
    """
    A data model for MIRI imaging apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - subarray: str[15]
        - eefraction: float32
        - radius: float32
        - apcorr: float32
        - skyin: float32
        - skyout: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirimg_apcorr.schema"


class MirLrsApcorrModel(ReferenceFileModel):
    """
    A data model for MIRI LRS apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - subarray: str[15]
        - wavelength: float32 1D array
        - nelem_wl: int16
        - size: uint8 1D array
        - nelem_size: int16
        - apcorr: float32 2D array
        - apcorr_err: float32 2D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirlrs_apcorr.schema"


class MirMrsApcorrModel(ReferenceFileModel):
    """
    A data model for MIRI MRS apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - wavelength: float32 1D array
        - radius: float32 2D array
        - apcorr: float32 2D array
        - apcorr_err: float32 2D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirmrs_apcorr.schema"


class NrcImgApcorrModel(ReferenceFileModel):
    """
    A data model for NIRCam imaging apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - eefraction: float32
        - radius: float32
        - apcorr: float32
        - skyin: float32
        - skyout: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrcimg_apcorr.schema"


class NrcWfssApcorrModel(ReferenceFileModel):
    """
    A data model for NIRCam WFSS apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - wavelength: float32 1D array
        - nelem_wl: int16
        - size: uint8 1D array
        - nelem_size: int16
        - apcorr: float32 2D array
        - apcorr_err: float32 2D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrcwfss_apcorr.schema"


class NisImgApcorrModel(ReferenceFileModel):
    """
    A data model for NIRISS imaging apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - eefraction: float32
        - radius: float32
        - apcorr: float32
        - skyin: float32
        - skyout: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nisimg_apcorr.schema"


class NisWfssApcorrModel(ReferenceFileModel):
    """
    A data model for NIRISS WFSS apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - wavelength: float32 1D array
        - nelem_wl: int16
        - size: uint8 1D array
        - nelem_size: int16
        - apcorr: float32 2D array
        - apcorr_err: float32 2D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/niswfss_apcorr.schema"


class NrsMosApcorrModel(ReferenceFileModel):
    """
    A data model for NIRSpec MOS apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - grating: str[15]
        - wavelength: float32 1D array
        - nelem_wl: int16
        - size: float32 2D array
        - nelem_size: int16
        - pixphase: float32 1D array
        - apcorr: float32 3D array
        - apcorr_err: float32 3D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrsmos_apcorr.schema"


class NrsIfuApcorrModel(ReferenceFileModel):
    """
    A data model for NIRSpec IFU apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - grating: str[15]
        - wavelength: float32 1D array
        - radius: float32 3D array
        - apcorr: float32 3D array
        - apcorr_err: float32 3D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrsifu_apcorr.schema"


class NrsFsApcorrModel(ReferenceFileModel):
    """
    A data model for NIRSpec Fixed-Slit apcorr reference files.

    Parameters
    __________
    apcorr_table : numpy table
        Aperture correction factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and aperture correction
        factors associated with those modes.

        - filter: str[12]
        - grating: str[15]
        - slit: str[15]
        - wavelength: float32 1D array
        - nelem_wl: int16
        - size: float32 2D array
        - nelem_size: int16
        - pixphase: float32 1D array
        - apcorr: float32 3D array
        - apcorr_err: float32 3D array

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrsfs_apcorr.schema"
