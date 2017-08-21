#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import CubeModel
from ..lib.catalog_utils import replace_suffix_ext
from .tso_photometry import tso_aperture_photometry


class TSOPhotometryStep(Step):
    """
    Perform circular aperture photometry on imaging Time Series
    Observations (TSO).

    Parameters
    -----------
    input : str or `CubeModel`
        A filename for either a FITS image or and association table or a
        `CubeModel`.
    """

    def process(self, input):
        with CubeModel(input) as model:
            # TODO:  need information about the actual source position in
            # TSO imaging mode (for all subarrays).
            # Meanwhile, this is a placeholder representing the geometric
            # center of the image.
            nint, ny, nx = model.data.shape
            xcenter = (ny - 1) / 2.
            ycenter = (ny - 1) / 2.

            # all radii are in pixel units
            if model.meta.instrument.pupil == 'WLP8':
                radius = 50
                radius_inner = 60
                radius_outer = 70
            else:
                radius = 3
                radius_inner = 4
                radius_outer = 5

            catalog = tso_aperture_photometry(model, xcenter, ycenter,
                                              radius, radius_inner,
                                              radius_outer)

            old_suffixes = ['calints', 'crfints']
            output_dir = self.search_attr('output_dir')
            cat_filepath = replace_suffix_ext(model.meta.filename,
                                              old_suffixes, 'phot',
                                              output_ext='ecsv',
                                              output_dir=output_dir)
            catalog.write(cat_filepath, format='ascii.ecsv', overwrite=True)
            self.log.info('Wrote TSO photometry catalog: {0}'.
                          format(cat_filepath))

        return catalog


if __name__ == '__main__':
    cmdline.step_script(tso_photometry_step)
