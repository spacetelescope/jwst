import warnings

from .multislit import MultiSlitModel


__all__ = ['MultiProductModel']


def MultiProductModel(*args, **kwargs):
    warnings.simplefilter('default')
    warnings.warn(message="MultiProductModel is deprecated and will be removed.  "
        "Use MultiSlitModel.", category=DeprecationWarning)
    return MultiSlitModel(*args, **kwargs)
