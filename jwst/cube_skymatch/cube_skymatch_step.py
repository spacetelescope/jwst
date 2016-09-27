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
        subtract = boolean(default=False) # subtract computed sky from cube data?

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

    def process(self, input):
        cube_models = datamodels.ModelContainer(input)

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

        # save background info in 'meta':
        for c in skymatch_cubes:
            self._set_bkg_meta(c.meta['original_cube_model'], c)

        return cube_models

    def _extend_schema(self, cube):
        bkg_schema = {
            "type": "object",
            "properties": {
                "meta": {
                    "type": "object",
                    "properties": {
                        "background": {
                            "type": "object",
                            "properties": {
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

        cube.extend_schema(bkg_schema)

    def _set_bkg_meta(self, cube_model, sky_cube):
        self._extend_schema(cube_model)

        cs_type = "image" if sky_cube.wcs is None else "world"
        coeffs = sky_cube.bkg_coeff.ravel().tolist()

        cube_model.meta.background.degree = sky_cube.bkg_degree
        cube_model.meta.background.refpoint = sky_cube.bkg_center
        cube_model.meta.background.cs_type = cs_type
        cube_model.meta.background.coefficients = coeffs


if __name__ == '__main__':
    cmdline.step_script(SkyMatchStep)



bkg_schema = {
    "type": "object",
    "properties": {
        "meta": {
            "type": "object",
            "properties": {
                "background": {
                    "type": "object",
                    "properties": {
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
                        "wcs": {"$ref": "http://stsci.edu/schemas/asdf/wcs/wcs-1.0.0"},
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
