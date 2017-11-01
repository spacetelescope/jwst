"""
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
__version__ = '0.7.1'
__vdate__ = '29-February-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
