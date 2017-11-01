""" This package provides support for image intensity subtraction and equalization
(matching) for MIRI images.

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)
import os
import logging


__docformat__ = 'restructuredtext'

from . mrs_imatch_step import MRSIMatchStep

__version__ = '0.7.1'
__vdate__ = '1-June-2017'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
