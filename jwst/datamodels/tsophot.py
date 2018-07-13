import warnings

from .reference import ReferenceFileModel
from .validate import ValidationWarning

__all__ = ['TsoPhotModel']

class TsoPhotModel(ReferenceFileModel):
    """
    A model for a reference file of type "tsophot".
    """
    schema_url = "tsophot.schema.yaml"
    reftype = "tsophot"

    def __init__(self, init=None, radii=None, **kwargs):
        super(TsoPhotModel, self).__init__(init=init, **kwargs)
        if radii is not None:
            self.radii = radii
        if init is None:
            self.populate_meta()

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

    def populate_meta(self):
        self.meta.instrument.name = "MIRI|NIRCAM|"
        self.meta.exposure.type = "MIR_IMAGE|NRC_TSIMAGE|"

    def to_fits(self):
        raise NotImplementedError("FITS format is not supported for this file.")

    def validate(self):
        super(TsoPhotModel, self).validate()
        try:
            assert len(self.radii) > 0
            assert self.meta.instrument.name in ["MIRI", "NIRCAM"]
            assert self.meta.exposure.type in ["MIR_IMAGE", "NRC_TSIMAGE"]
        except AssertionError as errmsg:
            if self._strict_validation:
                raise AssertionError(errmsg)
            else:
                warnings.warn(str(errmsg), ValidationWarning)
