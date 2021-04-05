from stcal.dynamicdq import dynamic_mask
from .dqflags import pixel
from .reference import ReferenceFileModel


__all__ = ['NirspecFlatModel', 'NirspecQuadFlatModel']


class NirspecFlatModel(ReferenceFileModel):
    """A data model for NIRSpec flat-field reference files.

    Parameters
    __________
    data : numpy float32 array
         NIRSpec flat-field reference data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error estimate

    wavelength : numpy table
         Table of wavelengths for image planes

    flat_table : numpy table
         Table for quickly varying component of flat field

    dq_def : numpy table
         DQ flag definitions
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nirspec_flat.schema"

    def __init__(self, init=None, **kwargs):
        super(NirspecFlatModel, self).__init__(init=init, **kwargs)

        if self.dq is not None or self.dq_def is not None:
            self.dq = dynamic_mask(self, pixel)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class NirspecQuadFlatModel(ReferenceFileModel):
    """A data model for NIRSpec flat-field files that differ by quadrant.

    Parameters
    __________
    quadrants.items.data : numpy float32 array


    quadrants.items.dq : numpy uint32 array


    quadrants.items.err : numpy float32 array


    quadrants.items.wavelength : numpy table
         Table of wavelengths for image planes

    quadrants.items.flat_table : numpy table
         Table for quickly varying component of flat field

    dq_def : numpy table
         DQ flag definitions
    """

    schema_url = "http://stsci.edu/schemas/jwst_datamodel/nirspec_quad_flat.schema"

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
