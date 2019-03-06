try:
    from .version import * # noqa: F403, F401
except ImportError:  # Not available for RTD
    __version_commit__ = 'unknown'
    __version__ = 'dev'
import sys

if sys.version_info < (3, 5):
    raise ImportError("JWST does not support Python 2.x, 3.0, 3.1, 3.2, 3.3 or 3.4."
                      "Beginning with JWST 0.9, Python 3.5 and above is required.")

from .transforms.jwextension import JWSTExtension
from .datamodels.extension import BaseExtension
