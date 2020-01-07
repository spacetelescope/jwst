from .reference import ReferenceFileModel

__all__ = ['FgsImgApcorrModel', 'MirImgApcorrModel',
           'NrcImgApcorrModel', 'NisImgApcorrModel']


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
