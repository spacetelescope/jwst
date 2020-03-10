""" skymatch

This package provides support for sky background subtraction and equalization
(matching).

"""
import logging

from .cube_skymatch_step import CubeSkyMatchStep

__author__ = 'Mihai Cara'


__all__ = ["CubeSkyMatchStep"]


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
