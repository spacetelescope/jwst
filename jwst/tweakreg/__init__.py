""" tweakreg

This package provides support for image alignment.

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)
from .tweakreg_step import TweakRegStep

import os
import logging

#from . import distortion

__docformat__ = 'restructuredtext'

from . import wcsutils
from . import imalign
from . import wcsimage
from . import matchutils
from . import tweakreg_step
from . import chelp
from . import linearfit
from . import simplewcs

__taskname__ = 'tweakreg'
__version__ = '0.7.1'
__vdate__ = '17-April-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def help():
    msg = \
"""
The tweakreg package contains the following tasks that allow users
perform WCS alignment.

tweakreg:
       align - primary task for performing image alignment
"""
    print(msg)
