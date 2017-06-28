#! /usr/bin/env python
# 
#  guider_cds.py

from __future__ import division
import time
import numpy as np
import logging

from .. import datamodels
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

BUFSIZE = 1024 * 30000  # 30Mb cache size for data section

def guider_cds(model, buffsize):
    """

    """

    # return new_model
