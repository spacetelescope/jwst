import warnings

from stdatamodels.validate import ValidationWarning
from stcal.dynamicdq import dynamic_mask

from .model_base import JwstDataModel
from .dqflags import pixel


__all__ = ['ReferenceFileModel']


class ReferenceFileModel(JwstDataModel):
    """
    A data model for reference tables

    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/referencefile.schema"

    def __init__(self, init=None, **kwargs):
        super(ReferenceFileModel, self).__init__(init=init, **kwargs)
        self._no_asdf_extension = True
        self.meta.telescope = "JWST"

    def validate(self):
        """
        Convenience function to be run when files are created.
        Checks that required reference file keywords are set.
        """
        to_fix = []
        to_check = ['description', 'reftype', 'author', 'pedigree', 'useafter']
        for field in to_check:
            if getattr(self.meta, field) is None:
                to_fix.append(field)
        if self.meta.instrument.name is None:
            to_fix.append('instrument.name')
        if self.meta.telescope != 'JWST':
            to_fix.append('telescope')
        if to_fix:
            self.print_err(f'Model.meta is missing values for {to_fix}')
        super().validate()

    def print_err(self, message):
        if self._strict_validation:
            raise ValueError(message)
        else:
            warnings.warn(message, ValidationWarning)


class ReferenceImageModel(ReferenceFileModel):
    """
    A data model for 2D reference images.

    Reference image data model.

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/referenceimage.schema"

    def __init__(self, init=None, **kwargs):
        super(ReferenceImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

        if self.hasattr('dq_def'):
            self.dq = dynamic_mask(self, pixel)


class ReferenceCubeModel(ReferenceFileModel):
    """
    A data model for 3D reference images

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/referencecube.schema"

    def __init__(self, init=None, **kwargs):
        super(ReferenceCubeModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class ReferenceQuadModel(ReferenceFileModel):
    """
    A data model for 4D reference images

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/referencequad.schema"

    def __init__(self, init=None, **kwargs):
        super(ReferenceQuadModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
