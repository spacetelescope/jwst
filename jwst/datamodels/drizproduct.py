import warnings

from .model_base import DataModel
from .image import ImageModel

__all__ = ['DrizProductModel']


class DrizProductModel(DataModel):
    """
    A data model for drizzle-generated products.

    Parameters
    __________
    data : numpy float32 array
         The science data

    con : numpy int32 array
         Drizzle Context array

    wht : numpy float32 array
         Drizzle Weight array
    """
    schema_url = "http://stsci.edu/schemas/jwst_datamodel/drizproduct.schema"

    @property
    def hdrtab(self):
        return self._extra_fits.HDRTAB.data

    @hdrtab.setter
    def hdrtab(self, v):
        self._extra_fits.HDRTAB.data = v

def DrizProductModel(*args, **kwargs):
    warnings.simplefilter('default')
    warnings.warn(message="DrizProduct is deprecated and will be removed.  "
        "Use ImageModel.", category=DeprecationWarning)
    return ImageModel(*args, **kwargs)
