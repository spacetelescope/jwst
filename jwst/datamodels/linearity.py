from stcal.dynamicdq import dynamic_mask
from .dqflags import pixel
from .reference import ReferenceFileModel


__all__ = ['LinearityModel']


class LinearityModel(ReferenceFileModel):
    """
    A data model for linearity correction information.

    Parameters
    __________
    coeffs : numpy float32 array
         Linearity coefficients

    dq : numpy uint32 array
         Data quality flags

    dq_def : numpy table
         DQ flag definitions
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/linearity.schema"

    def __init__(self, init=None, **kwargs):
        super(LinearityModel, self).__init__(init=init, **kwargs)

        self.dq = dynamic_mask(self, pixel)

        # Implicitly create arrays
        self.dq = self.dq

    def get_primary_array_name(self):
        """
        Returns the name "primary" array for this model, which
        controls the size of other arrays that are implicitly created.
        This is intended to be overridden in the subclasses if the
        primary array's name is not "data".
        """
        return 'coeffs'
