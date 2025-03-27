#!/usr/bin/env python
from stdatamodels.jwst.datamodels import CubeModel, GainModel
import numpy as np

from ..stpipe import Step
from ..lib.pipe_utils import is_tso
from .tso_photometry_centroiding import tso_photometry_centroiding
# from .tso_spectroscopy_centroiding import tso_spectroscopy_centroiding

__all__ = ['TSOCentroidingStep']


class TSOCentroidingStep(Step):
    """
    Estimate the centroid and PSF-width of a source in a TSO observation.

    Parameters
    ----------
    input : str or `CubeModel`
        Filename for a FITS image or association table, or a `CubeModel`.

    Returns
    -------
    catalog : tuple
        A tuple of three numpy arrays, the first containing the centroid
        results, the second containing the PSF width results, and the third
        containing the PSF flux results.
    """

    class_alias = "tso_centroiding"

    reference_file_types = ['gain']

    def process(self, input_data):

        # Open the input as a CubeModel
        with CubeModel(input_data) as model:

            # Get the gain reference file
            gain_filename = self.get_reference_file(model, 'gain')
            gain_model = GainModel(gain_filename)

            # Imaging
            if (model.meta.exposure.type == 'NRC_TSIMAGE' or
                    (model.meta.exposure.type == 'MIR_IMAGE' and is_tso(model))):
                catalog = tso_photometry_centroiding(model, gain_model)

            # Spectroscopy
            else:
                # catalog = tso_spectroscopy_centroiding(model, gain_model)

                # Just going to return zeros for now and change this later
                # FINDME: This won't work for MIRI which is rotated 90 degrees;
                #         will need to determine the dispersion axis from the
                #         header.
                catalog = (np.zeros((model.data.shape[0], model.data.shape[2])),
                           np.zeros((model.data.shape[0], model.data.shape[2])),
                           np.zeros((model.data.shape[0], model.data.shape[2])))

        return catalog
