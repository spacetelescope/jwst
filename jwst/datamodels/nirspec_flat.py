from .reference import ReferenceFileModel
from .dynamicdq import dynamic_mask

__all__ = ['NirspecFlatModel', 'NirspecQuadFlatModel']


class NirspecFlatModel(ReferenceFileModel):
    """A data model for NIRSpec flat-field reference files.

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data.  2-D or 3-D.

    dq: numpy array
        The data quality array.  2-D or 3-D.

    err: numpy array
        The error array.  2-D or 3-D.

    wavelength: numpy array
        The wavelength for each plane of the `data` array.  This will
        only be needed if `data` is 3-D.

    flat_table: numpy array
        A table of wavelengths and flat-field values, to specify the
        component of the flat field that can vary over a relatively short
        distance (can be pixel-to-pixel).
    """

    schema_url = "nirspec_flat.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(NirspecFlatModel, self).__init__(init=init, **kwargs)

        if self.dq is not None or self.dq_def is not None:
            self.dq = dynamic_mask(self)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class NirspecQuadFlatModel(ReferenceFileModel):
    """A data model for NIRSpec flat-field files that differ by quadrant.

    Parameters
    ----------
    init: any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.

    data: numpy array
        The science data.  2-D or 3-D.

    dq: numpy array
        The data quality array.  2-D or 3-D.

    err: numpy array
        The error array.  2-D or 3-D.

    wavelength: numpy array
        The wavelength for each plane of the `data` array.  This will
        only be needed if `data` is 3-D.

    flat_table: numpy array
        A table of wavelengths and flat-field values, to specify the
        component of the flat field that can vary over a relatively short
        distance (can be pixel-to-pixel).

    dq_def: numpy array
        The data quality definitions table.
    """

    schema_url = "nirspec_quad_flat.schema.yaml"

    def __init__(self, init=None, **kwargs):
        if isinstance(init, NirspecFlatModel):
            super(NirspecQuadFlatModel, self).__init__(init=None, **kwargs)
            self.update(init)
            self.quadrants.append(self.quadrants.item())
            self.quadrants[0].data = init.data
            self.quadrants[0].dq = init.dq
            self.quadrants[0].err = init.err
            self.quadrants[0].wavelength = init.wavelength
            self.quadrants[0].flat_table = init.flat_table
            self.quadrants[0].dq_def = init.dq_def
            self.quadrants[0].dq = dynamic_mask(self.quadrants[0])
            return

        super(NirspecQuadFlatModel, self).__init__(init=init, **kwargs)
