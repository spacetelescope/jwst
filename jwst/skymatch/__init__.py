""" skymatch

This package provides support for sky background subtraction and equalization
(matching).

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)
from .skymatch_step import SkyMatchStep

import os
import logging


__docformat__ = 'restructuredtext'

from . import skystatistics
from . import skymatch
from . import skyimage
from . import skymatch_step

__taskname__ = 'skymatch'
__version__ = '0.1.0'
__vdate__ = '29-February-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def help():
    msg = \
"""
The skymatch package contains the following tasks that allow users
perform sky level matching on user images.

skymatch:
       match - primary task for performing sky level matching on user images
"""
    print(msg)
