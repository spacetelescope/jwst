#! /usr/bin/env python
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..stpipe import Step, cmdline
from .. import datamodels
import logging
from .assign_wcs import load_wcs

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class AssignWcsStep(Step):
    """
    AssignWcsStep: Loads all WCS and distortion information for an exposure
    and stores it in the model meta data.
    """

    spec = """
    """

    """
    Reference file types:

    camera             NIRSPEC Camera model
    collimator         NIRSPEC Collimator Model
    disperser          Disperser parameters
    distortion         Spatial distortion model
    filteroffset       MIRI Imager fiter offsets
    fore               Transform through the NIRSPEC FORE optics
    fpa                Transform in the NIRSPEC FPA plane
    msa                Transformin the NIRSPEC MSA plane
    ote                Transform through the Optical Telescope Element
    specwcs            Wavelength calibration models
    regions            Stores location of the regions on the detector
    v2v3               Transform from MIRI instrument focal plane to V2V3 plane
    wavelengthrange    Typical wavelength ranges
    """
    reference_file_types = ['distortion', 'filteroffset', 'specwcs', 'regions',
                            'wavelengthrange', 'v2v3', 'camera', 'collimator',
                            'disperser', 'fore', 'fpa', 'msa', 'ote', 'ifupost',
                            'ifufore', 'ifuslicer']

    def process(self, input):
        reference_file_names = {}
        with datamodels.open(input) as input_model:
            # If input type is not supported, log warning, set to 'skipped', exit
            if not (isinstance(input_model, datamodels.ImageModel) or
                    isinstance(input_model, datamodels.CubeModel)):
                log.warning("Input dataset type is not supported, as it is")
                log.warning("neither ImageModel or CubeModel, so skipping step.")
                result = input_model.copy()
                result.meta.cal_step.assign_wcs = 'SKIPPED'
                return result

            for reftype in self.reference_file_types:
                reffile = self.get_reference_file(input_model, reftype)
                reference_file_names[reftype] = reffile if reffile else ""
            result = load_wcs(input_model, reference_file_names)

        return result


if __name__ == '__main__':
    cmdline.step_script(AssignWcsStep)
