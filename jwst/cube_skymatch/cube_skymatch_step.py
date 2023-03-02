#! /usr/bin/env python
"""
JWST pipeline step for sky matching.

:Authors: Mihai Cara

"""

import numpy as np

from ..stpipe import Step

from astropy.nddata.bitmask import (
    bitfield_to_boolean_mask,
    interpret_bit_flags,
)

from stdatamodels.jwst.datamodels.dqflags import pixel

from jwst.datamodels import ModelContainer

# LOCAL:
from .skymatch import match
from .skycube import SkyCube
from ..skymatch.skystatistics import SkyStats

__all__ = ['CubeSkyMatchStep']


class CubeSkyMatchStep(Step):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.

    """

    class_alias = "cube_skymatch"

    spec = """
        # General sky matching parameters:
        bkg_degree = integer(min=0, default=1) # Degree of the polynomial for background fitting
        subtract = boolean(default=False) # subtract computed sky from input1 cube data?
        subtract2d = boolean(default=True) # subtract computed sky from input2 2D data?

        # Sky statistics parameters:
        skystat = option('mode', 'median', 'mode', 'midpt', default='mode') # sky statistics
        dqbits = string(default=None) # "good" DQ bits
        lower = float(default=None) # Lower limit of "good" pixel values
        upper = float(default=None) # Upper limit of "good" pixel values
        nclip = integer(min=0, default=5) # number of sky clipping iterations
        lsigma = float(min=0.0, default=4.0) # Lower clipping limit, in sigma
        usigma = float(min=0.0, default=4.0) # Upper clipping limit, in sigma
        binwidth = float(min=0.0, default=0.1) # Bin width for 'mode' and 'midpt' `skystat`, in sigma
    """

    reference_file_types = []

    def process(self, input1, input2):
        cube_models = ModelContainer(input1)
        models2d = ModelContainer(input2)
        dqbits = interpret_bit_flags(self.dqbits, flag_name_map=pixel)

        # set sky statistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth
        )

        # At this moment running 'cube_skymatch' on images whose
        # background has been previously subtracted is not supported.
        # Raise an exception if background was subtracted:
        self._check_background(cube_models)
        self._check_background(models2d)
        self._reset_background(cube_models)
        self._reset_background(models2d)

        # create a list of SkyCubes:
        skycubes = []

        for cm in cube_models:
            # process weights and combine with DQ:
            if not hasattr(cm, 'weightmap') or cm.weightmap is None:
                weights = np.ones_like(cm.data, dtype=np.float64)
            else:
                weights = cm.weightmap.copy()

            if dqbits is not None:
                dq = bitfield_to_boolean_mask(
                    cm.dq,
                    self._dqbits,
                    good_mask_value=0,
                    dtype=bool
                )
                weights[dq] = 0.0

            wcs = cm.meta.wcs if hasattr(cm.meta, 'wcs') else None
            wcsinfo = cm.meta.wcsinfo if hasattr(cm.meta, 'wcsinfo') else None

            exptime = cm.meta.exposure.exposure_time
            if exptime is None:
                exptime = 1.0

            sc = SkyCube(
                data=cm.data,
                wcs=wcs,
                wcsinfo=wcsinfo,
                weights=weights,
                cube_weight=exptime,
                bkg_deg=self.bkg_degree,
                bkg_center=None,
                id=None,
                meta={'original_cube_model': cm}
            )

            skycubes.append(sc)

        skymatch_cubes, nsubspace = match(skycubes, subtract=self.subtract)

        if nsubspace > 1:
            self.log.warning("Not all cubes have been sky matched as "
                             "some of them do not overlap.")

        # save background info in 'meta' and subtract sky from 2D images
        # if requested:
        # model.meta.instrument.channel

        for c in skymatch_cubes:
            self._set_cube_bkg_meta(c.meta['original_cube_model'], c)
            model2d, channel = _find_associated_2d_image(
                c.meta['original_cube_model'], models2d
            )

            if model2d is None:
                continue

            self._set_model2d_bkg_meta(
                c.meta['original_cube_model'],
                model2d,
                channel
            )

            if self.subtract2d:
                self._apply_sky_2d(model2d, channel)

        return cube_models, models2d

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

    def _set_cube_bkg_meta(self, model, sky_cube):
        np = len(model.meta.background.polynomial_info)
        if np > 1 or (np == 1 and model.meta.background.polynomial_info[0].channel):
            raise ValueError("Cube's 'polynomial_info' must be empty or "
                             "contain at most one element with a 'channel' "
                             "property set to None.")

        pinfo = {
            'degree': list(sky_cube.bkg_degree),
            'refpoint': list(sky_cube.bkg_center),
            'cs_type': "image" if sky_cube.wcs is None else "world",
            'coefficients': sky_cube.bkg_coeff.ravel().tolist(),
            'wcs': sky_cube.wcs
        }

        model.meta.background.subtract = self.subtract
        model.meta.background.level = None
        if np == 0:
            model.meta.background.polynomial_info.append(pinfo)
        else:
            model.meta.background.polynomial_info[0] = pinfo

    def _set_model2d_bkg_meta(self, model3d, model2d, channel):
        channel = str(channel)

        # see if meta for this channel exist and
        # find the index of this channel:
        index = _find_channel_bkg_index(model2d, channel)

        cubebkg = model3d.meta.background.polynomial_info[0]
        bkgmeta = {
            "degree": list(cubebkg.degree),
            "refpoint": list(cubebkg.refpoint),
            "cs_type": cubebkg.cs_type,
            "wcs": cubebkg.wcs,
            "coefficients": list(cubebkg.coefficients),
            "channel": channel
        }

        model2d.meta.background.subtract = self.subtract2d
        model2d.meta.background.level = None
        if index is None:
            model2d.meta.background.polynomial_info.append(bkgmeta)
        else:
            model2d.meta.background.polynomial_info[index] = bkgmeta

    def _apply_sky_2d(self, model2d, channel):
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
        x, y = self._get_2d_pixgrid(model2d, channel)
        x = x.ravel()
        y = y.ravel()

        # convert to RA/DEC:
        ra, dec, lam = model2d.meta.wcs(x.astype(dtype=float),
                                        y.astype(dtype=float))

        # some pixels may be NaNs and so throw them out:
        m = np.logical_and(
            np.logical_and(np.isfinite(ra), np.isfinite(dec)),
            np.isfinite(lam)
        )
        x = x[m]
        y = y[m]
        ra = ra[m]
        dec = dec[m]
        lam = lam[m]

        if cs_type == "image":
            raise ValueError("Polynomials must be defined in world coords in "
                             "order to be able to perform sky subtraction "
                             "on 2D DataModel.")

        elif cs_type != "world":
            raise ValueError("Unsupported background polynomial's cs_type'.")

        # compute background values:
        ra -= refpt[0]
        dec -= refpt[1]
        lam -= refpt[2]
        bkg = np.polynomial.polynomial.polyval3d(ra, dec, lam, c)

        # subtract background:
        model2d.data[y, x] -= bkg

    def _get_2d_pixgrid(self, model2d, channel):
        # TODO: the code in this function is experimental and most likely will
        # will need revising at a later time. At this moment, I was told we
        # cannot use WCS domain to find the range of pixel indices that
        # belong to a given channel. Therefore, for now we have this hardcoded
        # in this function.
        y, x = np.indices((1024, 512))

        if channel in ['1', '3']:
            return x + 4, y
        else:
            return x + 516, y


def _find_associated_2d_image(cube, models2d):
    # TODO: the code in this function is experimental and most likely will
    # not work. I created it mostly as a placeholder so that the rest of the
    # code could be developed. However, at this moment it is not clear how
    # a 2D image from which a cube was created could be identified.
    """
    Returns a tuple consisting of a 2D image and channel number from which
    a cube was created.

    """
    fn = cube.meta.filename
    ch = int(cube.meta.instrument.channel)

    for m in models2d:
        if m.meta.filename == fn:
            return m, ch

    return None, None


def _find_channel_bkg_index(model2d, channel):
    """
    Return the index of the background subschema corresponding to a given
    channel.

    """
    channel = str(channel)
    index = None
    for k, m in enumerate(model2d.meta.background.polynomial_info):
        if m.channel == channel:
            index = k
    return index
