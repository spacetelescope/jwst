from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['CubeModel']


class CubeModel(model_base.DataModel):
    """
    A data model for 3D image cubes.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The science data.  3-D.

    dq : numpy array
        The data quality array.  3-D.

    err : numpy array
        The error array.  3-D
    """
    schema_url = "cube.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, zeroframe=None,
                 relsens=None, area=None, **kwargs):
        super(CubeModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        if zeroframe is not None:
            self.zeroframe = zeroframe

        if relsens is not None:
            self.relsens = relsens

        if area is not None:
            self.area = area

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

