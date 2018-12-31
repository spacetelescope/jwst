import warnings
import sys
import traceback

from .model_base import DataModel
from .dynamicdq import dynamic_mask
from .validate import ValidationWarning

__all__ = ['ReferenceFileModel']


class ReferenceFileModel(DataModel):
    """
    A data model for reference tables

    """
    schema_url = "referencefile.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(ReferenceFileModel, self).__init__(init=init, **kwargs)
        self._no_asdf_extension = True
        self.meta.telescope = "JWST"

    def validate(self):
        """
        Convenience function to be run when files are created.
        Checks that required reference file keywords are set.
        """
        try:
            assert self.meta.description is not None
            assert self.meta.telescope == 'JWST'
            assert self.meta.reftype is not None
            assert self.meta.author is not None
            assert self.meta.pedigree is not None
            assert self.meta.useafter is not None
            assert self.meta.instrument.name is not None
        except AssertionError:
            if self._strict_validation:
                raise
            else:
                tb = sys.exc_info()[-1]
                tb_info = traceback.extract_tb(tb)
                text = tb_info[-1][-1]
                warnings.warn(text, ValidationWarning)

        super(ReferenceFileModel, self).validate()


class ReferenceImageModel(ReferenceFileModel):
    """
    A data model for 2D reference images

Reference image data model

    Parameters
    __________
    data : numpy float32 array
         The science data

    dq : numpy uint32 array
         Data quality array

    err : numpy float32 array
         Error array
    """
    schema_url = "referenceimage.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(ReferenceImageModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

        if self.hasattr('dq_def'):
            self.dq = dynamic_mask(self)


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
    schema_url = "referencecube.schema.yaml"

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
    schema_url = "referencequad.schema.yaml"

    def __init__(self, init=None, **kwargs):
        super(ReferenceQuadModel, self).__init__(init=init, **kwargs)

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
