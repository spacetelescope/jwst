import warnings
from .model_base import DataModel
from .validate import ValidationWarning

__all__ = ['ReferenceFileModel']


class ReferenceFileModel(DataModel):
    """
    A data model for reference tables

    Parameters
    ----------
    init : any
        Any of the initializers supported by `~jwst.datamodels.DataModel`.
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
            assert (self.meta.telescope == 'JWST')
            assert self.meta.reftype is not None
            assert self.meta.author is not None
            assert self.meta.pedigree is not None
            assert self.meta.useafter is not None
            assert self.meta.instrument.name is not None
        except AssertionError as errmsg:
            if self._strict_validation:
                raise AssertionError(errmsg)
            else:
                warnings.warn(str(errmsg), ValidationWarning)

        super(ReferenceFileModel, self).validate()


class ReferenceImageModel(ReferenceFileModel):
    """
    A data model for 2D reference images

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
    """
    schema_url = "referenceimage.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, **kwargs):
        super(ReferenceImageModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err


class ReferenceCubeModel(ReferenceFileModel):
    """
    A data model for 3D reference images

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
    """
    schema_url = "referencecube.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, **kwargs):
        super(ReferenceCubeModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err

class ReferenceQuadModel(ReferenceFileModel):
    """
    A data model for 4D reference images

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
    """
    schema_url = "referencequad.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, **kwargs):
        super(ReferenceQuadModel, self).__init__(init=init, **kwargs)

        if data is not None:
            self.data = data

        if dq is not None:
            self.dq = dq

        if err is not None:
            self.err = err

        # Implicitly create arrays
        self.dq = self.dq
        self.err = self.err
