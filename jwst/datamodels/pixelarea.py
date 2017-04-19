from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['PixelAreaModel', 'NirspecAreaModel']


class PixelAreaModel(model_base.DataModel):
    """
    A data model for the pixel area map
    """
    schema_url = "pixelarea.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        super(PixelAreaModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

class NirspecAreaModel(model_base.DataModel):
    """
    A data model for the NIRSpec pixel area reference file

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    slit_table : numpy array
        A table-like object containing row selection criteria made up
        of the slit id and the pixel area values associated with
        the slits.

        slit_id: str[15]
        pixarea: float32

    mos_table : numpy array
        A table-like object containing row selection criteria made up
        of MOS shutter parameters and the pixel area values associated
        with the shutters.

        quadrant: int16
        shutter_x: int16
        shutter_y: int16
        pixarea: float32

    ifu_table : numpy array
        A table-like object containing row selection criteria made up
        of IFU slice id and the pixel area values associated
        with the slices.

        slice_id: int16
        pixarea: float32
    """
    schema_url = "nirspec_area.schema.yaml"

    def __init__(self, init=None, slit_table=None, mos_table=None,
        ifu_table=None, **kwargs):
        super(NirspecAreaModel, self).__init__(init=init, **kwargs)

        if slit_table is not None:
            self.slit_table = slit_table

        if mos_table is not None:
            self.mos_table = mos_table

        if ifu_table is not None:
            self.ifu_table = ifu_table

