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

    ntables = len(input.spec)
    fluxsums = []
    times = []

    # The input should contain one row per integration for each spectral
    # order.  NIRISS SOSS data can contain up to three orders.
    norders = 0                 # number of different spectral orders
    sporders = []               # list of spectral order numbers
    ntables_order = []          # number of tables for each spectral order
    prev_spectral_order = -999
    for i in range(ntables):
        # The following assumes that all rows for a given spectral order
        # are contiguous.
        spectral_order = input.spec[i].spectral_order
        if spectral_order is None:
            norders = 1
            sporders = [0]
            ntables_order = [ntables]
            break
        if spectral_order != prev_spectral_order:
            sporders.append(spectral_order)
            prev_spectral_order = spectral_order
            ntables_order.append(1)
            norders = len(sporders)
        else:
            ntables_order[norders - 1] += 1

    nints = max(ntables_order)
    log.debug("norders = %d, sporders = %s, ntables_order = %s",
              norders, str(sporders), str(ntables_order))

    # Compute the flux sum for each integration in the input
    for i in range(ntables):
        fluxsums.append(input.spec[i].spec_table['FLUX'].sum())

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
    dt_arr = np.zeros(ntables, dtype=np.float64)
    dt = (input.meta.exposure.group_time * (input.meta.exposure.ngroups + 1))
    j0 = 0
    for k in range(norders):
        ntables_current = ntables_order[k]
        j1 = j0 + ntables_current
        dt_arr[j0 : j1] = np.arange(1, 1 + ntables_current) * dt - (dt / 2.)
        j0 += ntables_current
    int_dt = TimeDelta(dt_arr, format='sec')

    # Compute the absolute time at the mid-point of each integration
    # Note that this won't be correct if the input has been segmented and
    # the current file is not the first segment.
    int_times = (Time(input.meta.exposure.start_time, format='mjd') + int_dt)

    # Store the times and flux sums in the table
    tbl['MJD'] = int_times.mjd
    tbl['whitelight_flux'] = fluxsums

    return tbl
