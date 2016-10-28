from __future__ import absolute_import

from os.path import (
    abspath,
    dirname,
    join
)
__version__ = '0.7.0-beta.1'


# Utility
def libpath(filepath):
    '''Return the full path to the module library.'''

    return join(dirname(abspath(__file__)),
                'lib',
                filepath)

from .association import *
from .exceptions import *
from .generate import *
from .pool import *
from .registry import *
