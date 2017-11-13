""" tweakreg

This package provides support for image alignment.

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

import os
import logging

__docformat__ = 'restructuredtext'

__taskname__ = 'tweakreg'
__version__ = '0.8.0'
__vdate__ = '10-November-2017'
__author__ = 'Mihai Cara'

from .tweakreg_step import TweakRegStep
from . import wcsutils
from . import imalign
from . import wcsimage
from . import matchutils
from . import tweakreg_step
from . import chelp
from . import linearfit
from . import tpcorr


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
