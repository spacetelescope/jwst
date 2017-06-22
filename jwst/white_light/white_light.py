from __future__ import absolute_import

import time
import logging

import numpy as np
from collections import OrderedDict
from astropy.table import QTable
from astropy.time import Time, TimeDelta
from ..datamodels import dqflags

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def white_light(input):

    nints = input.meta.exposure.nints
    fluxsums = []
    times = []

    # Compute the flux sum for each integration in the input
    for i in np.arange(nints):
        fluxsums.append(input.spec[i].spec_table['flux'].sum())

    # Populate meta data for the output table
    tbl_meta = OrderedDict()
    tbl_meta['instrument'] = input.meta.instrument.name
    tbl_meta['detector'] = input.meta.instrument.detector
    tbl_meta['exp_type'] = input.meta.exposure.type
    tbl_meta['subarray'] = input.meta.subarray.name
    tbl_meta['filter'] = input.meta.instrument.filter
    tbl_meta['pupil'] = input.meta.instrument.pupil
    tbl_meta['target_name'] = input.meta.target.catalog_name
    tbl_meta['number_of_integrations'] = nints

    # Create the output table
    tbl = QTable(meta=tbl_meta)

    # Compute the delta time of each integration
    dt = (input.meta.exposure.group_time * (input.meta.exposure.ngroups + 1))
    dt_arr = (np.arange(1, 1 + input.meta.exposure.nints) * dt - (dt / 2.))
    int_dt = TimeDelta(dt_arr, format='sec')

    # Compute the absolute time at the mid-point of each integration
    int_times = (Time(input.meta.exposure.start_time, format='mjd') + int_dt)

    # Store the times and flux sums in the table
    tbl['MJD'] = int_times.mjd
    tbl['whitelight_flux'] = fluxsums

    return tbl
