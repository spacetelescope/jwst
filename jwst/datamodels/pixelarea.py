from .reference import ReferenceFileModel


__all__ = ['PixelAreaModel', 'NirspecSlitAreaModel', 'NirspecMosAreaModel',
           'NirspecIfuAreaModel']


class PixelAreaModel(ReferenceFileModel):
    """
    A data model for the pixel area map

    Parameters
    __________
    data : numpy float32 array
         The pixel area array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/pixelarea.schema"



class NirspecSlitAreaModel(ReferenceFileModel):
    """
    A data model for the NIRSpec fixed-slit pixel area reference file

    Parameters
    __________
    area_table : numpy table
         NIRSpec fixed-slit pixel area table
         A table-like object containing row selection criteria made up
         of the slit id and the pixel area values associated with
         the slits.

         - slit_id: str[15]
         - pixarea: float32
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nirspec_area_slit.schema"


class NirspecMosAreaModel(ReferenceFileModel):
    """
    A data model for the NIRSpec MOS pixel area reference file

    Parameters
    __________
    area_table : numpy table
         NIRSpec MOS pixel area table
         A table-like object containing row selection criteria made up
         of MOS shutter parameters and the pixel area values associated
         with the shutters.

         - quadrant: int16
         - shutter_x: int16
         - shutter_y: int16
         - pixarea: float32
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nirspec_area_mos.schema"


class NirspecIfuAreaModel(ReferenceFileModel):
    """
    A data model for the NIRSpec IFU pixel area reference file

    Parameters
    __________
    area_table : numpy table
         NIRSpec IFU pixel area table
         A table-like object containing row selection criteria made up
         of IFU slice id and the pixel area values associated
         with the slices.

         - slice_id: int16
         - pixarea: float32
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nirspec_area_ifu.schema"
