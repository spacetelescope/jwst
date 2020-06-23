"""V1 Calculation based on time and engineering database info
"""
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

import jwst.datamodels as dm
import jwst.lib.set_telescope_pointing as stp

__all__ = ['v1_calculate_from_models']


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

    # Calculate V1 for all sources.
    for source in sources:
        with dm.open(source) as model:
            time_start = model.meta.exposure.start_time
            time_end = model.meta.exposure.end_time

            _, vinfo = stp.calc_wcs_from_time(time_start, time_end, **calc_wcs_from_time_kwargs)

            v1_dict['time_start'].append(time_start)
            v1_dict['time_end'].append(time_end)
            v1_dict['vinfo'].append(vinfo)

    # Format and return.
    v1_table = Table(v1_dict)
    return v1_table
