from stcal.dynamicdq import dynamic_mask
from .dqflags import pixel
from .reference import ReferenceFileModel


__all__ = ['FgsImgPhotomModel', 'MirImgPhotomModel', 'MirLrsPhotomModel',
           'MirMrsPhotomModel', 'NrcImgPhotomModel', 'NrcWfssPhotomModel',
           'NisImgPhotomModel', 'NisSossPhotomModel', 'NisWfssPhotomModel',
           'NrsFsPhotomModel', 'NrsMosPhotomModel']


class FgsImgPhotomModel(ReferenceFileModel):
    """
    A data model for FGS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - photmjsr: float32
        - uncertainty: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/fgsimg_photom.schema"


class MirImgPhotomModel(ReferenceFileModel):
    """
    A data model for MIRI imaging photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

       - filter: str[12]
       - subarray: str[15]
       - photmjsr: float32
       - uncertainty: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirimg_photom.schema"


class MirLrsPhotomModel(ReferenceFileModel):
    """
    A data model for MIRI LRS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

       - filter: str[12]
       - subarray: str[15]
       - photmjsr: float32
       - uncertainty: float32
       - nelem: int16
       - wavelength: float32[*]
       - relresponse: float32[*]
       - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirlrs_photom.schema"


class MirMrsPhotomModel(ReferenceFileModel):
    """
    A data model for MIRI MRS photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        An array-like object containing the pixel-by-pixel conversion values
        in units of (MJy / pixel) / (DN / sec).

    err : numpy array
        An array-like object containing the uncertainties in the conversion
        values, in the same units as the data array.

    dq : numpy array
        An array-like object containing bit-encoded data quality flags,
        indicating problem conditions for values in the data array.

    dq_def : numpy array
        A table-like object containing the data quality definitions table.

    pixsiz : numpy array
        An array-like object containing pixel-by-pixel size values, in units of
        square arcseconds (arcsec^2).
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/mirmrs_photom.schema"

    def __init__(self, init=None, **kwargs):
        super(MirMrsPhotomModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self, pixel)


class NrcImgPhotomModel(ReferenceFileModel):
    """
    A data model for NIRCam imaging photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[12]
        - photmjsr: float32
        - uncertainty: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrcimg_photom.schema"


class NrcWfssPhotomModel(ReferenceFileModel):
    """
    A data model for NIRCam WFSS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - order: int16
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[*]
        - relresponse: float32[*]
        - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrcwfss_photom.schema"


class NisImgPhotomModel(ReferenceFileModel):
    """
    A data model for NIRISS imaging photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[12]
        - photmjsr: float32
        - uncertainty: float32

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nisimg_photom.schema"


class NisWfssPhotomModel(ReferenceFileModel):
    """
    A data model for NIRISS WFSS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - order: int16
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[*]
        - relresponse: float32[*]
        - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/niswfss_photom.schema"


class NisSossPhotomModel(ReferenceFileModel):
    """
    A data model for NIRISS SOSS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[15]
        - order: int16
        - photmj: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[*]
        - relresponse: float32[*]
        - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nissoss_photom.schema"


class NrsFsPhotomModel(ReferenceFileModel):
    """
    A data model for NIRSpec Fixed-Slit photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - grating: str[15]
        - slit: str[15]
        - photmj: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[*]
        - relresponse: float32[*]
        - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrsfs_photom.schema"


class NrsMosPhotomModel(ReferenceFileModel):
    """
    A data model for NIRSpec MOS and IFU photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - grating: str[15]
        - photmj: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[*]
        - relresponse: float32[*]
        - reluncertainty: float32[*]

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nrsmos_photom.schema"
