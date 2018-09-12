"""Setup default and environment """

# Take version from the upstream package
from .. import __version__


# Utility
def libpath(filepath):
    '''Return the full path to the module library.'''
    from os.path import (
        abspath,
        dirname,
        join
    )
    return join(dirname(abspath(__file__)),
                'lib',
                filepath)

from .association import *
from .exceptions import *
from .association_io import *
from .generate import *
from .pool import *
from .registry import *
from .load_asn import load_asn
