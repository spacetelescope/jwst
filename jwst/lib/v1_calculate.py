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
    siaf = stp.SIAF(v2_ref=0., v3_ref = 0., v3yangle = 0., vparity = 1.)

    # Calculate V1 for all sources.
    for source in sources:
        with dm.open(source) as model:
            obsstart = model.meta.exposure.start_time
            obsend = model.meta.exposure.end_time

            _, vinfos = stp.calc_wcs_from_time(obsstart, obsend, siaf=siaf, **calc_wcs_from_time_kwargs)

    # Format and return.
    v1_table = Table(v1_dict)
    return v1_table
