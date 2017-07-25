#! /usr/bin/env python
"""
JWST pipeline step for image intensity matching for MIRI images.

:Authors: Mihai Cara

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

import numpy as np

from .. stpipe import Step, cmdline
from .. import datamodels
from .. wiimatch.match import *
from .. cube_build import CubeBuildStep

try:
    from stsci.tools.bitmask import bitfield_to_boolean_mask
    from stsci.tools.bitmask import interpret_bit_flags
except ImportError:
    from stsci.tools.bitmask import bitmask2mask as bitfield_to_boolean_mask
    from stsci.tools.bitmask import interpret_bits_value as interpret_bit_flags


__all__ = ['MRSIMatchStep']


class MRSIMatchStep(Step):
    """
    MRSIMatchStep: Subtraction or equalization of sky background in MIRI MRS science images.

    """

    spec = """
        # General sky matching parameters:
        bkg_degree = integer(min=0, default=1) # Degree of the polynomial for background fitting
        subtract = boolean(default=False) # subtract computed sky from 'images' cube data?

    """

    reference_file_types = []

    def process(self, images):
        all_models2d = datamodels.ModelContainer(images, persist=True)

        chm = {}

        for m in all_models2d:
            ch = m.meta.instrument.channel
            if ch not in chm:
                chm[ch] = datamodels.ModelContainer()
            chm[ch].append(m)

        # check that channel combinations are reasonable, in particular that
        # a given channel does not appear in multiple channel lists.
        # For example, input image list cannot contain some images with
        # channel='12' and other images with channel='13' (channel 1 combines
        # with *both* 2 & 3).
        single_ch = {}
        for ch in chm.keys():
            for c in map(int, list(ch)):
                if c in single_ch:
                    raise ValueError("Channel '{:d}' appears in multiple "
                                     "combinations of channels.".format(c))
                single_ch[c] = chm[ch]

        # At this moment running 'cube_skymatch' on images whose
        # background has been previously subtracted is not supported.
        # Raise an exception if background was subtracted:
        for m2d in chm.values():
            self._check_background(m2d)

        # reset background in a separate loop after all checks passed:
        for m2d in chm.values():
            self._reset_background(m2d)

        # check that 'degree' has the same length as the number of dimensions
        # in image arrays:
        if hasattr(self.bkg_degree, '__iter__'):
            if len(self.bkg_degree) != 3:
                raise ValueError("The length of 'bkg_degree' parameter must "
                                 "be 3 or 'bkg_degree' must be an integer.")
            degree = tuple([int(d) for d in self.bkg_degree])

        else:
            intdeg = int(self.bkg_degree)
            degree = (intdeg, intdeg, intdeg)

        # match background for images from a single channel
        for c in sorted(single_ch.keys()):
            matched_models = _match_models(single_ch[c], channel=str(c),
                                           degree=degree)
            if self.subtract:
                for m in matched_models:
                    _apply_sky_2d(m, channel=str(c))

        return images

    def _check_background(self, models):
        # see if 'cube_skymatch' step was previously run and raise
        # an exception as 'cube_skymatch' cannot be run twice on the
        # same data:
        for m in models:
            if m.meta.background.subtracted is None:
                if m.meta.background.level is not None:
                    # report inconsistency:
                    raise ValueError("Background level was set but the "
                                     "'subtracted' property is undefined "
                                     "(None).")

            elif m.meta.background.subtracted:
                raise ValueError("'cube_skymatch' step cannot be run on "
                                 "data whose background has been previously "
                                 "subtracted.")

    def _reset_background(self, models):
        for m in models:
            del m.meta.background


def _apply_sky_2d(model2d, channel):
    """ Apply/subtract sky from 2D image data. """

    index = _find_channel_bkg_index(model2d, channel)
    if index is None:
        raise ValueError("Background data for channel '{}' not present in "
                         "2D model '{}'"
                         .format(channel, model2d.meta.filename))

    # get may parameters of the background polynomial:
    bkgmeta = model2d.meta.background.polynomial_info[index]
    degree = tuple(bkgmeta.degree)
    degree_p1 = tuple((i + 1 for i in degree))
    c = np.reshape(list(bkgmeta.coefficients), degree_p1)
    refpt = tuple(bkgmeta.refpoint)

    cs_type = bkgmeta.cs_type

    # get pixel grid for sky computations:
    x, y = _get_2d_pixgrid(model2d, channel)
    x = x.ravel()
    y = y.ravel()

    # convert to RA/DEC:
    r, d, l = model2d.meta.wcs(x.astype(dtype=np.float),
                               y.astype(dtype=np.float))

    # some pixels may be NaNs and so throw them out:
    m = np.logical_and(
        np.logical_and(np.isfinite(r), np.isfinite(d)),
        np.isfinite(l)
    )
    x = x[m]
    y = y[m]
    r = r[m]
    d = d[m]
    l = l[m]

    if cs_type == "image":
        raise ValueError("Polynomials must be defined in world CS in "
                         "order to be able to perform sky subtraction "
                         "on 2D DataModel.")

    elif cs_type != "world":
        raise ValueError("Unsupported background polynomial's 'cs_type'.")

    # compute background values:
    r -= refpt[0]
    d -= refpt[1]
    l -= refpt[2]
    bkg = np.polynomial.polynomial.polyval3d(r, d, l, c)

    # subtract background:
    model2d.data[y, x] -= bkg


def _get_2d_pixgrid(model2d, channel):
    # TODO: the code in this function is experimental and most likely will
    # will need revising at a later time. At this moment, I was told we
    # cannot use WCS domain to find the range of pixel indices that
    # belong to a given channel. Therefore, for now we have this hardcoded
    # in this function.
    y, x = np.indices((1024, 512))

    if channel in ['1', '3']:
        return (x + 4, y)
    else:
        return (x + 516, y)


def _match_models(models, channel, degree, center=None, center_cs='image'):
    # create a list of cubes:
    cbs = CubeBuildStep()
    cbs.channel = str(channel)
    cbs.band = 'ALL'
    cbs.single = True
    cube_models = cbs.process(models)
    if len(cube_models) != len(models):
        raise RuntimeError("The number of generated cube models does not "
                           "match the number of input 2D images.")

    # retrieve WCS (all cubes must have identical WCS so we use the first):
    meta = cube_models[0].meta

    if hasattr(meta, 'wcs'):
        wcs = meta.wcs
    else:
        raise ValueError("Cubes build from input 2D images do not contain WCS."
                         " Unable to proceed.")

    wcsinfo = meta.wcsinfo if hasattr(meta, 'wcsinfo') else None
    if wcsinfo is not None and (wcsinfo.crval1 is None or
                                wcsinfo.crval2 is None or
                                wcsinfo.crval3 is None):
        raise ValueError("'wcsinfo' cannot have its 'crvaln' set to None.")

    # set center of the coordinate system to CRVAL if available:
    if center is None and wcsinfo is not None:
        center = (wcsinfo.crval1, wcsinfo.crval2, wcsinfo.crval3)
        center_cs = 'world'

    # build lists of data, masks, and sigmas (weights)
    image_data = []
    mask_data = []
    sigma_data = []

    for cm in cube_models:
        #TODO: at this time it is not clear that data should be weighted
        #      by exptime the way it is done below and possibly should be
        #      revised later.
        exptime = cm.meta.exposure.exposure_time
        if exptime is None:
            exptime = 1.0

        # process weights and create masks:
        if not hasattr(cm, 'weightmap') or cm.weightmap is None:
            weights = np.ones_like(cm.data, dtype=np.float64)
            sigmas = weights / np.sqrt(exptime)
            mask = np.ones_like(weights, dtype=np.uint8)
            mask_data.append(mask)

        else:
            weights = cm.weightmap.copy()
            eps = np.finfo(weights.dtype).tiny
            bad_data = weights < eps
            weights[bad_data] = eps # in order to avoid runtime warnings
            sigmas = 1.0 / np.sqrt(exptime * weights)
            mask = np.logical_not(bad_data).astype(np.uint8)
            mask_data.append(mask)

        image_data.append(cm.data)
        sigma_data.append(sigmas)

    bkg_poly_coef, mat, _, _, effc, cs = match_lsq(
        images=image_data,
        masks=mask_data,
        sigmas=sigma_data,
        degree=degree,
        center=center,
        image2world=wcs.__call__,
        center_cs=center_cs,
        ext_return=True
    )

    if cs != 'world':
        raise RuntimeError("Unexpected coordinate system.")

    #TODO: try to identify if all images overlap
    #if nsubspace > 1:
        #self.log.warning("Not all cubes have been sky matched as "
                         #"some of them do not overlap.")

    # save background info in 'meta' and subtract sky from 2D images
    # if requested:
    ##### model.meta.instrument.channel

    # set 2D models' background meta info:
    for im, poly in zip(models, bkg_poly_coef):
        im.meta.background.subtracted = False
        im.meta.background.polynomial_info.append(
            {
                'degree': degree,
                'refpoint': center,
                'cs_type': 'world',
                'wcs': wcs,
                'coefficients': poly.ravel().tolist(),
                'channel': channel
            }
        )

    return models


def _find_channel_bkg_index(model2d, channel):
    """
    Return the index of the background subschema corrsponding to a given
    channel.

    """
    channel = str(channel)
    index = None
    for k, m in enumerate(model2d.meta.background.polynomial_info):
        if m.channel == channel:
            index = k
    return index


if __name__ == '__main__':
    cmdline.step_script(SkyMatchStep)
