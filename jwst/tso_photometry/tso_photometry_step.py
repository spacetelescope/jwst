#!/usr/bin/env python

from ..stpipe import Step, cmdline
from ..datamodels import CubeModel
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
        with CubeModel(input) as datamodel:
            # TODO:  need information about the actual source position in
            # TSO imaging mode (for all subarrays).
            # Meanwhile, this is a placeholder representing the geometric
            # center of the image.
            nint, ny, nx = datamodel.data.shape
            xcenter = (ny - 1) / 2.
            ycenter = (ny - 1) / 2.

            # all radii are in pixel units
            if datamodel.meta.instrument.pupil == 'WLP8':
                radius = 50
                radius_inner = 60
                radius_outer = 70
            else:
                radius = 3
                radius_inner = 4
                radius_outer = 5

            catalog = tso_aperture_photometry(datamodel, xcenter, ycenter,
                                              radius, radius_inner,
                                              radius_outer)
            catalog_format = 'ecsv'
            fmt = 'ascii.ecsv'
            catalog_filename = datamodel.meta.filename.replace(
                '.fits', '_cat.{0}'.format(catalog_format))

            catalog.write(catalog_filename, format=fmt)
            self.log.info('Wrote source catalog: {0}'.
                          format(catalog_filename))

        return


if __name__ == '__main__':
    cmdline.step_script(tso_photometry_step)
