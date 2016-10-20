#! /usr/bin/env python
"""
JWST pipeline step for sky matching.

:Authors: Mihai Cara

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""
from __future__ import (absolute_import, division, unicode_literals,
                        print_function)

import numpy as np

from ..stpipe import Step, cmdline
from .. import datamodels

from stsci.tools import bitmask

#LOCAL:
from . skymatch import match
from . skycube import SkyCube
from ..skymatch.skystatistics import SkyStats

__all__ = ['CubeSkyMatchStep']


class CubeSkyMatchStep(Step):
    """
    SkyMatchStep: Subtraction or equalization of sky background in science images.

    """

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
        cube_models = datamodels.ModelContainer(input1)
        models2d = datamodels.ModelContainer(input2)

        dqbits = bitmask.interpret_bits_value(self.dqbits)

        # set sky stattistics:
        self._skystat = SkyStats(
            skystat=self.skystat,
            lower=self.lower,
            upper=self.upper,
            nclip=self.nclip,
            lsig=self.lsigma,
            usig=self.usigma,
            binwidth=self.binwidth
        )

        # create a list of SkyCubes:
        skycubes = []

        for cm in cube_models:

            # process weights and combine with DQ:
            if not hasattr(cm, 'weightmap') or cm.weightmap is None:
                weights = np.ones_like(cm.data, dtype=np.float64)
            else:
                weights = cm.weightmap.copy()

            if dqbits is not None:
                dq = bitmask.bitmask2mask(bitmask=cm.dq,
                                          ignore_bits=self._dqbits,
                                          good_mask_value=0, dtype=np.bool)
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
        ##### model.meta.instrument.channel

        for c in skymatch_cubes:
            self._set_cube_bkg_meta(c.meta['original_cube_model'], c)
            model2d, channel = _find_associated_2d_image(
                c.meta['original_cube_model'], models2d
            )

            if model2d is None:
                continue

            self._set_model2d_bkg_meta(c, model2d, channel)

            if self.subtract2d:
                self._apply_sky_2d(model2d, channel)

        return cube_models

    def _set_cube_bkg_meta(self, model, sky_cube):
        _extend_cube_schema(model)

        cs_type = "image" if sky_cube.wcs is None else "world"
        coeffs = sky_cube.bkg_coeff.ravel().tolist()

        model.meta.background.degree = sky_cube.bkg_degree
        model.meta.background.refpoint = sky_cube.bkg_center
        model.meta.background.cs_type = cs_type
        model.meta.background.coefficients = coeffs

    def _set_model2d_bkg_meta(self, model3d, model2d, channel):
        _extend_2d_schema(model2d)
        channel = str(channel)

        # see if meta for this channel exist and
        # find the index of this channel:
        index = _find_channel_bkg_index(model2d, channel)

        cubebkg = model3d.meta.background
        bkgmeta = {
            "degree": cubebkg.degree,
            "refpoint": cubebkg.refpoint,
            "cs_type": cubebkg.cs_type,
            "wcs": cubebkg.wcs,
            "coefficients": cubebkg.coefficients,
            "channel": channel
        }

        if index is None:
            model2d.meta.background.append(bkgmeta)
        else:
            model2d.meta.background[index] = bkgmeta

    def _apply_sky_2d(self, model2d, channel):
        """ Apply/subtract sky from 2D image data. """

        index = _find_channel_bkg_index(model2d, channel)
        if index is None:
            raise ValueError("Background data for channel '{}' not present in "
                             "2D model '{}'"
                             .format(channel, model2d.meta.filename))

        # get may parameters of the background polynomial:
        bkgmeta = model2d.meta.background[index]
        degree = bkgmeta.degree
        degree_p1 = tuple((i + 1 for i in degree))
        c = np.reshape(bkgmeta.coefficients, degree_p1)
        refpt = bkgmeta.refpoint

        cs_type = bkgmeta.cs_type

        # get pixel grid for sky computations:
        x, y = self._get_2d_pixgrid(model2d, channel)
        x = x.ravel()
        y = y.ravel()

        # convert to RA/DEC:
        r, d, l = model2d.met.wcs(x.astype(dtype=np.float),
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
            raise ValueError("Unsupported background polynomial's cs_type'.")

        # compute background values:
        r -= refpt[0]
        d -= refpt[1]
        l -= refpt[2]
        bkg = np.polynomial.polynomial.polyval3d(r, d, l, c)

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
            return (x + 4, y)
        else:
            return (x + 516, y)


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
            return (m, ch)

    return None, None


def _find_channel_bkg_index(model2d, channel):
    """
    Return the index of the background subschema corrsponding to a given
    channel.

    """
    channel = str(channel)
    index = None
    for k, m in enumerate(model2d.meta.background):
        if m.channel == channel:
            index = k
    return index


def _extend_cube_schema(self, model):
    bkg_schema = {
        "type": "object",
        "properties": {
            "meta": {
                "type": "object",
                "properties": {
                    "background": {
                        "type": "object",
                        "properties": {
                            "channel": {
                                "title": "Channel for which background is "
                                         "applicable",
                                "type": "string"
                            },

                            "degree": {
                                "title": "Max degree in each coordinate",
                                "type": "array",
                                "items": [
                                    {
                                      "type": "number"
                                    }
                                ]
                            },

                            "refpoint": {
                                "title": "Reference coordinates",
                                "type": "array",
                                "items": [
                                    {
                                      "type": "number"
                                    }
                                ]
                            },

                            "cs_type": {
                                "title": "Type of coordinates.",
                                "type": "string",
                                "enum": ["world", "image"]
                            },

                            "wcs": {"$ref": "http://stsci.edu/schemas/asdf"
                                    "/wcs/wcs-1.0.0"},

                            "coefficients": {
                                "title": "Polynomial Coefficients",
                                "type": "array",
                                "items": [
                                    {
                                      "type": "number"
                                    }
                                ]
                            }
                        }
                    }
                }
            }
        }
    }

    model.extend_schema(bkg_schema)


def _extend_2d_schema(model):
    bkg_schema = {
        "type": "object",
        "properties": {
            "meta": {
                "type": "object",
                "properties": {
                    "background": {
                        "title": "List of background objects for each channel",
                        "type": "array",
                        "items": [
                            {
                              "type": "object",
                              "properties": {

                                  "degree": {
                                      "title": "Max degree in each "
                                               "coordinate",
                                      "type": "array",
                                      "items": [
                                          {
                                            "type": "number"
                                          }
                                      ]
                                  },

                                  "refpoint": {
                                      "title": "Reference coordinates",
                                      "type": "array",
                                      "items": [
                                          {
                                            "type": "number"
                                          }
                                      ]
                                  },

                                  "cs_type": {
                                      "title": "Type of coordinates.",
                                      "type": "string",
                                      "enum": ["world", "image"]
                                  },

                                  "wcs": {"$ref": "http://stsci.edu/"
                                          "schemas/asdf/wcs/wcs-1.0.0"},

                                  "coefficients": {
                                      "title": "Polynomial Coefficients",
                                      "type": "array",
                                      "items": [
                                          {
                                            "type": "number"
                                          }
                                      ]
                                  },

                                  "channel": {
                                      "title": "Channel ID",
                                      "type": "string"
                                  }
                              }
                            }
                        ]
                    }
                }
            }
        }
    }

    model.extend_schema(bkg_schema)


if __name__ == '__main__':
    cmdline.step_script(SkyMatchStep)
