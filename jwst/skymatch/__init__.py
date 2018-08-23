"""
This package provides support for sky background subtraction and equalization
(matching).

"""
import logging
from .skymatch_step import SkyMatchStep

__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["SkyMatchStep"]
