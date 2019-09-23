from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['PhotomModel']


class PhotomModel(ReferenceFileModel):
    """
    A base class for photometric reference file models.

    Parameters
    __________
    phot_table : numpy table
         Photometric flux conversion factors table
    """
    schema_url = "photom.schema"


class NircamPhotomModel(PhotomModel):
    """
    A data model for NIRCam photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[12]
        - order: int16
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[3000]
        - relresponse: float32[3000]

    """
    schema_url = "nircam_photom.schema"


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
        - order: int16
        - photmjsr: float32
        - uncertainty: float32

    """
    schema_url = "nrcimg_photom.schema"


class NirissPhotomModel(PhotomModel):
    """
    A data model for NIRISS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - pupil: str[12]
        - order: int16
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[5000]
        - relresponse: float32[5000]

    """
    schema_url = "niriss_photom.schema"


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
    schema_url = "nisimg_photom.schema"


class NirspecPhotomModel(PhotomModel):
    """
    A data model for NIRSpec imaging, IFU, and MOS photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - grating: str[12]
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[150]
        - relresponse: float32[150]
        - reluncertainty: float32[150]

    """
    schema_url = "nirspec_photom.schema"


class NirspecFSPhotomModel(PhotomModel):
    """
    A data model for NIRSpec Fixed-Slit (FS) photom reference files.

    Parameters
    __________
    phot_table : numpy table
        Photometric flux conversion factors table
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - filter: str[12]
        - grating: str[12]
        - slit: str[12]
        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[150]
        - relresponse: float32[150]
        - reluncertainty: float32[150]

    """
    schema_url = "nirspecfs_photom.schema"

    def __init__(self, init=None, **kwargs):
        super(NirspecFSPhotomModel, self).__init__(init=init, **kwargs)


class MiriImgPhotomModel(PhotomModel):
    """
    A data model for MIRI imager photom reference files.

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
       - wavelength: float32[500]
       - relresponse: float32[500]
       - relresperror: float32[500]

    """
    schema_url = "miriimg_photom.schema"


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
    schema_url = "mirimg_photom.schema"


class MiriMrsPhotomModel(PhotomModel):
    """
    A data model for MIRI MRS photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        An array-like object containing the pixel-by-pixel conversion values
        in units of DN / sec / mJy / pixel.

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
    schema_url = "mirmrs_photom.schema"

    def __init__(self, init=None, **kwargs):
        super(MiriMrsPhotomModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self)


class FgsPhotomModel(PhotomModel):
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
        - nelem: int16
        - wavelength: float32[5000]
        - relresponse: float32[5000]

    """
    schema_url = "fgs_photom.schema"


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
    schema_url = "fgsimg_photom.schema"
