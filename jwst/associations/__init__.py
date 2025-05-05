"""
Association generator and rule definitions.

The association generator takes a list of items and an Association Pool, and it
creates sub-lists of those items depending on each item's attributes. The
association rules define how the sub-lists are created.
"""

# Take version from the upstream package
from jwst import __version__


# Utility
def libpath():
    """
    Return the full path to the module library.

    Returns
    -------
    Path
        Path to the module library.
    """
    from pathlib import Path

    return Path(__file__).parent / "lib"


from .association import *
from .association_io import *
from .exceptions import *
from .generator import *
from .lib.process_list import *
from .pool import *
from .registry import *
from .load_asn import load_asn
from .main import *
