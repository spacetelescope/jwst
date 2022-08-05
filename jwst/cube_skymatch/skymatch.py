"""
A module that provides functions for matching sky in overlapping cubes.

:Authors: Mihai Cara

"""

import logging
import copy
from datetime import datetime


__all__ = ['match']


__author__ = 'Mihai Cara'


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def match(skycubes, subtract=False):
    """
    A function to compute and/or "equalize" sky background in input cubes.

    .. note::
        Sky matching ("equalization") is possible only for **overlapping**
        cubes.

    .. note::
        This function may modify input `SkyCube` objects. Make a copy of the
        input objects if their preservation is important.


    Parameters
    ----------
    skycubes : list of SkyCubes
        A list of of :py:class:`SkyCube` objects whose sky needs to be
        matched.

    subtract : bool, optional
        Indicates whether to subtract computed sky differences from cube data.


    Returns
    -------

    out_skycubes : list of SkyCube
        A list of `SkyCube` each containing sky background polynomial
        coefficients and, if requested, background-subtracted data.

        .. note::
            Output list is not in the same order as input list but rather it
            follows the order in which cubes have been matched with the first
            cube being a "reference" cube.

    nsubspace : int
        Number of subspaces. A value larger than one indicates that not all
        images overlap.


    """
    function_name = match.__name__

    # Time it
    runtime_begin = datetime.now()

    # deep-copy input so that input is not affected by changes:
    # skycubes = copy.deepcopy(skycubes)
    skycubes = [s for s in skycubes]

    nsubspace = 0
    out_skycubes = []

    if len(skycubes) < 2:
        # too few images: nothing to match to
        log.info("Need at least two images for sky background matching.")

        # there are no more overlapping cubes. Add the remaining cubes to
        # the output and return:
        out_skycubes += skycubes
        nsubspace += len(skycubes)
        skycubes = []

    while len(skycubes) > 0:
        # find a pair of cubes with largest overlap and pick one of them as the
        # "reference" cube that we will "grow" to form a reference mosaic:

        # find largest overlap
        c1, c2, om = max_overlap_pair(skycubes)
        log.info("Overlap measure: {}".format(om))

        if c1 is None:
            if nsubspace == 0:
                log.warning("Input images do not overlap. Cannot match sky "
                            "background.")
            break

        out_skycubes.append(c1)
        nsubspace += 1

        # make a deep copy of c1 without deep-copying its 'meta':
        c1_meta = c1.meta
        c1.meta = {}

        # next 4 lines are a workaround to this bug:
        # https://github.com/STScI-JWST/jwst/issues/284 in DataModel
        # TODO: no need to pull out and later re-insert wcs and wcsinfo
        # once the bug is fixed.
        c1_wcs = c1._wcs
        c1_wcsinfo = c1._wcsinfo
        c1._wcs = None
        c1._wcsinfo = None

        mosaic = copy.deepcopy(c1)
        c1.meta = c1_meta

        # continuation of the workaround to DataModel bug.
        # TODO: remove this once the bug is fixed.
        c1._wcs = c1_wcs
        c1._wcsinfo = c1_wcsinfo

        # fit background
        c2.fit_background(mosaic)
        mosaic.combine_with_other(c2)
        out_skycubes.append(c2)

        # find cubes that overlap with the mosaic & match them
        while len(skycubes) > 0:
            c, om = max_overlap_cube(mosaic, skycubes)
            log.info("Overlap measure: {}".format(om))
            if c is None:
                # mosaic does not overlap with any of the remaining cubes
                break

            c.fit_background(mosaic)
            mosaic.combine_with_other(c)
            out_skycubes.append(c)

    # there are no more overlapping cubes. Add the remaining cubes to
    # the output and return:
    out_skycubes += skycubes
    nsubspace += len(skycubes)

    # subtract computed background from cube data if requested:
    if subtract:
        for c in out_skycubes:
            c.subtract_sky()

    # log running time:
    runtime_end = datetime.now()
    log.info(" ")
    log.info("***** {:s}.{:s}() ended on {}"
             .format(__name__, function_name, runtime_end))
    log.info("***** {:s}.{:s}() TOTAL RUN TIME: {}"
             .format(__name__, function_name, runtime_end - runtime_begin))
    log.info(" ")

    return out_skycubes, nsubspace


def max_overlap_pair(cubes):
    """
    Find, in the input list of cubes, the pair that has the largest overlap
    and remove these two cubes from the input list.

    """
    n = len(cubes)
    if n == 1:
        return None, None, 0.0  # leave single cube in the input list

    im = -1
    jm = -1
    om = -1.0  # max overlap
    for i in range(n):
        for j in range(i + 1, n):
            o = cubes[i].overlap_measure(cubes[j])
            if o > om:
                om = o
                im = i
                jm = j

    if om <= 0.0:
        return None, None, 0.0  # leave non-overlapping cubes in the input list

    if im < jm:
        jm -= 1
    c1 = cubes.pop(im)
    c2 = cubes.pop(jm)

    return c1, c2, om


def max_overlap_cube(refcube, cubes):
    n = len(cubes)

    if n < 1:
        return None, 0.0

    im = -1
    om = -1.0

    for i in range(n):
        o = refcube.overlap_measure(cubes[i])
        if o > om:
            om = o
            im = i

    if om > 0.0:
        c = cubes.pop(im)
    else:
        c = None  # leave non-overlapping cubes in the input list

    return c, om
