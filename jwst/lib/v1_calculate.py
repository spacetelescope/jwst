"""V1 Calculation based on time and engineering database info
"""
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from astropy.table import Table

import jwst.datamodels as dm
import jwst.lib.set_telescope_pointing as stp

__all__ = ['v1_calculate_from_models', 'v1_calculate_over_time']


def v1_calculate_from_models(sources, **calc_wcs_from_time_kwargs):
    """Calculate V1 over the time period for the given models

    Parameters
    ----------
    sources : [File-like or jwst.datamodels.Datamodel[...]]

    calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to `calc_wcs_from_time`

    Returns
    -------
    v1_table : astropy.table.Table
        Table of V1 pointing
    """
    # Initialize structures.
    v1_dict = defaultdict(list)
    siaf = stp.SIAF(v2_ref=0., v3_ref=0., v3yangle=0., vparity=1.)

    # Calculate V1 for all sources.
    for source in sources:
        with dm.open(source) as model:
            obsstart = model.meta.exposure.start_time
            obsend = model.meta.exposure.end_time

            obstimes, _, vinfos = stp.calc_wcs_over_time(
                obsstart, obsend, siaf=siaf, **calc_wcs_from_time_kwargs
            )
        sources = [source]*len(obstimes)
        v1_dict['source'] += sources
        v1_dict['obstime'] += obstimes
        v1_dict['v1'] += vinfos

    # Format and return.
    v1_table = Table(v1_dict)
    return v1_table


def v1_calculate_over_time(obsstart, obsend, **calc_wcs_from_time_kwargs):
    """Calculate V1 over the given time period

    Parameters
    ----------
    obsstart, obbsend : float
        The MJD start and end time to search for pointings.

    calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to `calc_wcs_from_time`

    Returns
    -------
    v1_table : astropy.table.Table
        Table of V1 pointing
    """
    # Initialize structures.
    siaf = stp.SIAF(v2_ref=0., v3_ref=0., v3yangle=0., vparity=1.)

    # Calculate V1 for all sources.
    obstimes, _, vinfos = stp.calc_wcs_over_time(
        obsstart, obsend, siaf=siaf, **calc_wcs_from_time_kwargs
    )
    v1_dict = dict()
    v1_dict['source'] = ['time range']*len(obstimes)
    v1_dict['obstime'] = obstimes
    v1_dict['v1'] = vinfos

    # Format and return.
    v1_table = Table(v1_dict)
    return v1_table
