"""Association Generator

The Association Generator takes a list of items, an Association Pool, and
creates sub-lists of those items depending on each item's attributes. How the
sub-lists are created is defined by Association Rules.

For more, see the :ref:`documentation overview <asn-overview>`.

"""

# Take version from the upstream package
from .. import __version__ # noqa F401


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

from .association import * # noqa F402,F403
from .association_io import * # noqa F402,F403
from .exceptions import * # noqa F402,F403
from .generator import * # noqa F402,F403
from .lib.process_list import * # noqa F402,F403
from .pool import * # noqa F402,F403
from .registry import * # noqa F402,F403
from .load_asn import load_asn # noqa F402,F403
from .main import * # noqa F402,F403
