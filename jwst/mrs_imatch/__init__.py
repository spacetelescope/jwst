"""Support for image intensity subtraction and equalization (matching) for MIRI images."""

import logging

from .mrs_imatch_step import MRSIMatchStep

__all__ = ["MRSIMatchStep"]

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
