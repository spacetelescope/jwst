""" skymatch

This package provides support for sky background subtraction and equalization
(matching).

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)
import os
import logging


__docformat__ = 'restructuredtext'

from . import skymatch
from . import skycube
from . import cube_skymatch_step

__taskname__ = 'cube_skymatch'
__version__ = '0.1.0'
__vdate__ = '21-Sept-2016'
__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def help():
    msg = \
"""
The skymatch package contains the following tasks that allow users
perform sky level matching on user images.

skymatch:
       match - primary task for performing sky level matching on user images.
       apply_match - subtract sky based on computed polynomials stored in meta.
"""
    print(msg)
