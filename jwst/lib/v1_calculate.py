"""V1 Calculation based on time and engineering database info
"""
from collections import defaultdict
import logging

from astropy.table import Table

from . import set_telescope_pointing as stp
from . import siafdb
from .. import datamodels as dm

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = ['v1_calculate_from_models', 'v1_calculate_over_time', 'simplify_table']


def v1_calculate_from_models(sources, siaf_path=None, **calc_wcs_from_time_kwargs):
    """Calculate V1 over the time period for the given models

    Parameters
    ----------
    sources : [File-like or jwst.datamodels.Datamodel[...]]
        The datamodels to get timings other header parameters from.

    siaf_path : None or file-like
        The path to the SIAF database. If none, the default used
        by the `pysiaf` package is used. See `SiafDb` for more
        information.

    calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to `calc_wcs_from_time`

    Returns
    -------
    v1_table : astropy.table.Table
        Table of V1 pointing
    """
    # Initialize structures.
    v1_dict = defaultdict(list)
    siaf = siafdb.SIAF(v2_ref=0., v3_ref=0., v3yangle=0., vparity=1.)
    t_pars = stp.TransformParameters(siaf=siaf, **calc_wcs_from_time_kwargs)

    # Calculate V1 for all sources.
    with siafdb.SiafDb(siaf_path) as siaf_db:
        t_pars.siaf_db = siaf_db
        for source in sources:
            with dm.open(source) as model:
                obsstart = model.meta.exposure.start_time
                obsend = model.meta.exposure.end_time
                guide_star_wcs = stp.WCSRef(
                    model.meta.guidestar.gs_ra,
                    model.meta.guidestar.gs_dec,
                    model.meta.guidestar.gs_v3_pa_science
                )

                t_pars.guide_star_wcs = guide_star_wcs
                obstimes, _, vinfos = stp.calc_wcs_over_time(obsstart, obsend, t_pars)

            sources = [source] * len(obstimes)
            v1_dict['source'] += sources
            v1_dict['obstime'] += obstimes
            v1_dict['v1'] += vinfos

    # Format and return.
    v1_table = Table(v1_dict, meta=t_pars.as_reprdict())
    return v1_table


def v1_calculate_over_time(obsstart, obsend, siaf_path=None, **calc_wcs_from_time_kwargs):
    """Calculate V1 over the given time period

    Parameters
    ----------
    obsstart, obsend : float
        The MJD start and end time to search for pointings.

    siaf_path : None or file-like
        The path to the SIAF database. If none, the default used
        by the `pysiaf` package is used. See `SiafDb` for more
        information.

    calc_wcs_from_time_kwargs : dict
        Keyword arguments to pass to `calc_wcs_from_time`

    Returns
    -------
    v1_table : astropy.table.Table
        Table of V1 pointing
    """
    # Initialize structures.
    siaf = siafdb.SIAF(v2_ref=0., v3_ref=0., v3yangle=0., vparity=1.)
    t_pars = stp.TransformParameters(siaf=siaf, **calc_wcs_from_time_kwargs)

    # Calculate V1 for all sources.
    with siafdb.SiafDb(siaf_path) as siaf_db:
        t_pars.siaf_db = siaf_db
        obstimes, _, vinfos = stp.calc_wcs_over_time(obsstart, obsend, t_pars)
    v1_dict = dict()
    v1_dict['source'] = ['time range'] * len(obstimes)
    v1_dict['obstime'] = obstimes
    v1_dict['v1'] = vinfos

    # Format and return.
    v1_table = Table(v1_dict, meta=t_pars.as_reprdict())
    return v1_table


def simplify_table(v1_table):
    """Convert pure object-based table to ASCII/Human-friendly

    The tables as produced by the `v1_calculate` functions use native objects.
    For instance, the "obstime" column contains `astopy.time.Time` objects and
    "v1" is the `jwst.lib.set_telescope_pointing.WCSREF` object

    This routine converts such objects to strings or Python-native builtin objects.

    Parameters
    ----------
    v1_table: astropy.table.Table
        V1 table as produced by `v1_calculate` functions.

    Returns
    -------
    formatted: astropy.table.Table
        Reformatted table.
    """
    source_formatted = [str(v) for v in v1_table['source']]
    obstime_formatted = v1_table['obstime'].isot
    ras, decs, pa_v3s = list(map(list, zip(*v1_table['v1'])))

    formatted = Table(
        [source_formatted, obstime_formatted, ras, decs, pa_v3s],
        names=('source', 'obstime', 'ra', 'dec', 'pa_v3'),
        meta=v1_table.meta
    )
    return formatted
