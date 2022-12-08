#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from .profile_fitting import pixel_replace_profile_fit
__all__ = ["PixelReplaceStep"]


class PixelReplaceStep(Step):
    """
    PixelReplaceStep: Module for replacing flagged bad pixels with an estimate
    of their flux, prior to spectral extraction.

    Attributes
    ----------
    n_adjacent_cols : int
        Number of adjacent columns (on either side of column containing a bad pixel) to use in
        creation of source profile, in cross-dispersion direction. The total number of columns
        used in the profile will be twice this number; on array edges, take adjacent columns until
        this number is reached.
    """

    class_alias = "pixel_replace"

    spec = """
        n_adjacent_cols = integer(default=3)    # Number of adjacent columns to use in creation of profile
    """

    def process(self, input):
        """Execute the step.

        Parameters
        ----------
        input: JWST data model

        Returns
        -------
        JWST data model
            This will be `input_model` if the step was skipped; otherwise,
            it will be a model containing data arrays with estimated fluxes
            for any bad pixels, now flagged as TO-BE-DETERMINED (DQ bit 7?).
        """

        # If more than one method ends up being utilized, may be convenient to
        # have a selector

        with datamodels.open(input) as input_model:

            result = input_model.copy()
            # If more than one 2d spectrum exists in input, call replacement
            # for each spectrum
            if isinstance(input_model, datamodels.MultiSlitModel):
                self.log.debug('Input is a MultiSlitModel.')

            # MIRI_LRS-FIXEDSLIT comes in ImageModel
            elif isinstance(input_model, datamodels.ImageModel):
                # MIRI LRS mode - any others?
                self.log.debug('Input is a ImageModel.')
                result =
            '''
            elif isinstance(input_model, datamodels.CubeModel):
                # It's a 3-D multi-integration model
                self.log.debug('Input is a CubeModel for a multiple integ. file')
            elif isinstance(input_model, datamodels.IFUCubeModel):
                self.log.debug('Input is an IFUCubeModel')
            elif isinstance(input_model, datamodels.SlitModel):
                # NRS_BRIGHTOBJ and MIRI LRS fixed-slit (resampled) modes
                self.log.debug('Input is a SlitModel')
            '''
            else:
                self.log.error(f'Input is a {str(type(input_model))}, ')
                self.log.error('which was not expected for pixel_replace')
                raise Exception
                #self.log.error('pixel_replace will be skipped.')
                #input_model.meta.cal_step.pixel_replace = 'SKIPPED'
                #return input_model

            self.record_step_status(result, 'master_background', success=False)

