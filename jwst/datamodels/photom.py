from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['PhotomModel']


class PhotomModel(ReferenceFileModel):
    """
    A base class for photometric reference file models.
    """
    schema_url = "photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(PhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


class NircamPhotomModel(PhotomModel):
    """
    A data model for NIRCam photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
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
    schema_url = "nircam_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(NircamPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


class NirissPhotomModel(PhotomModel):
    """
    A data model for NIRISS photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
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
    schema_url = "niriss_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(NirissPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


class NirspecPhotomModel(PhotomModel):
    """
    A data model for NIRSpec imaging, IFU, and MOS photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
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
    schema_url = "nirspec_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(NirspecPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


class NirspecFSPhotomModel(PhotomModel):
    """
    A data model for NIRSpec Fixed-Slit (FS) photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
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
    schema_url = "nirspecfs_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(NirspecFSPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


class MiriImgPhotomModel(PhotomModel):
    """
    A data model for MIRI imaging photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
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

    """
    schema_url = "mirimg_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(MiriImgPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table


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
    schema_url = "mirmrs_photom.schema.yaml"

    def __init__(self, init=None, data=None, err=None, dq=None, dq_def=None,
                 pixsiz=None, **kwargs):
        super(MiriMrsPhotomModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if err is not None:
            self.err = err

        if dq is not None:
            self.dq = dq

        if dq_def is not None:
            self.dq_def = dq_def

        if pixsiz is not None:
            self.pixsiz = pixsiz

        self.dq = dynamic_mask(self)

class FgsPhotomModel(PhotomModel):
    """
    A data model for FGS photom reference files.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    phot_table : numpy array
        A table-like object containing row selection criteria made up
        of instrument mode parameters and photometric conversion
        factors associated with those modes.

        - photmjsr: float32
        - uncertainty: float32
        - nelem: int16
        - wavelength: float32[5000]
        - relresponse: float32[5000]

    """
    schema_url = "fgs_photom.schema.yaml"

    def __init__(self, init=None, phot_table=None, **kwargs):
        super(FgsPhotomModel, self).__init__(init=init, **kwargs)

        if phot_table is not None:
            self.phot_table = phot_table
