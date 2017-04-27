from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['ImageModel']


class ImageModel(model_base.DataModel):
    """
    A data model for 2D images.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data : numpy array
        The science data.

    dq : numpy array
        The data quality array.

    err : numpy array
        The error array.

    relsens : numpy array
        The relative sensitivity table.

    relsens2d: numpy array
        The relative sensitivty 2D array.
    """
    schema_url = "image.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, relsens=None,
                 relsens2d=None, zeroframe=None, area=None, **kwargs):
        super(ImageModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if relsens is not None:
            self.relsens = relsens

        if relsens2d is not None:
            self.relsens2d = relsens2d

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if area is not None:
            self.area = area

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
