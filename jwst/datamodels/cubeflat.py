from __future__ import absolute_import, unicode_literals, division, print_function

from . import model_base


__all__ = ['CubeFlatModel']


class CubeFlatModel(model_base.DataModel):
    """
    A data model for 3D image cubes for flat-field reference data.

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst_lib.models.DataModel`.

    data : numpy array
        The array containing flat-field values.  3-D.

    wavelength : numpy array
        The array containing wavelength values.  3-D.

    dq : numpy array
        The data quality array.  3-D.

    dq_def : numpy array
        The data quality definitions table.
    """
    schema_url = "cubeflat.schema.yaml"

    def __init__(self, init=None, data=None, wavelength=None, dq=None,
                 dq_def=None, **kwargs):
        super(CubeFlatModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if wavelength is not None:
            self.wavelength = wavelength

        if dq is not None:
            self.dq = dq

        if dq_def is not None:
            self.dq_def = dq_def

        # Implicitly create array
        self.dq = self.dq
