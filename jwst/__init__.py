try:
    from .version import * # noqa: F403, F401
except ImportError:  # Not available for RTD
    __version_commit__ = 'unknown'
    __version__ = 'dev'
import sys

if sys.version_info < (3, 5):
    raise ImportError("JWST requires Python 3.5 and above.")
