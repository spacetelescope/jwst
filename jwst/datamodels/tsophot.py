import warnings
import sys
import traceback

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

    def on_save(self, path=None):
        self.meta.reftype = self.reftype

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
                raise
            else:
                tb = sys.exc_info()[-1]
                tb_info = traceback.extract_tb(tb)
                text = tb_info[-1][-1]
                warnings.warn(text, ValidationWarning)
