""" This package provides support for image intensity subtraction and equalization
(matching) for MIRI images.

"""
import os
import logging


__docformat__ = 'restructuredtext'

from . mrs_imatch_step import MRSIMatchStep

__version__ = '0.9.3'
__vdate__ = '1-June-2017'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
