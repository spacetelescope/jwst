"""Association Generator

The Association Generator takes a list of items, an Association Pool, and
creates sub-lists of those items depending on each item's attributes. How the
sub-lists are created is defined by Association Rules.

For more, see the :ref:`documentation overview <asn-overview>`.

"""

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
from .association_io import *
from .exceptions import *
from .generate import *
from .lib.process_list import *
from .pool import *
from .registry import *
from .load_asn import load_asn
from .main import *
