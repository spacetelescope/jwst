"""Process JWST data with Python."""

import re
from importlib.metadata import version

__all__ = ["__version__", "__version_commit__"]

__version__ = version(__name__)

if _match := re.match(r".*\+g(\w+)", __version__):
    __version_commit__ = _match.groups()[0]
    # clean up namespace
    del _match
else:
    __version_commit__ = ""
