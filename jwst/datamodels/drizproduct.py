import warnings

from .image import ImageModel

__all__ = ['DrizProductModel']


def DrizProductModel(*args, **kwargs):
    warnings.simplefilter('default')
    warnings.warn(message="DrizProduct is deprecated and will be removed.  "
                  "Use ImageModel.", category=DeprecationWarning)
    return ImageModel(*args, **kwargs)
