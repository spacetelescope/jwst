"""
This package provides support for sky background subtraction and equalization
(matching).

"""
import logging

from .skymatch_step import SkyMatchStep
from . import skystatistics
from . import skymatch
from . import skyimage
from . import skymatch_step


__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

__all__ = ["skyimage", "skymatch", "skymatch_step", "skystatistics",
           "SkyMatchStep"]
