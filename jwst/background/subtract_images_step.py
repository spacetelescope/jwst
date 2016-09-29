#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import subtract_images

class SubtractImagesStep(Step):
    """
    SubtractImagesStep:  Subtract two exposures from one
    another to accomplish background subtraction.
    """

    spec = """
    """

    def process(self, input1, input2):

        """
        Short Summary
        -------------
        Subtract the background signal from a JWST data model by
        subtracting a background image from it.

        Parameters
        ----------
        input1: JWST data model
            input science data model to be background-subtracted

        input2: JWST data model
            background data model

        Returns
        -------
        result: JWST data model
            background-subtracted science data model
        """

        # First, determine what kind of input model has been provided
        model1 = datamodels.open(input1)

        if isinstance(model1, datamodels.CubeModel):
            self.log.debug('Input is a CubeModel')
        elif isinstance(model1, datamodels.ImageModel):
            self.log.debug('Input is an ImageModel')
        elif isinstance(model1, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')

        # Assume that the second input model is always Image or MultiSlit with
        # a single image, which is safe to open as MultiSlit for either case
        #model2 = datamodels.MultiSlitModel(input2)
        model2 = datamodels.ImageModel(input2)

        # Call the subtraction routine
        result = subtract_images.subtract(model1, model2)

        # Set the step status indicator in the output model
        result.meta.cal_step.back_sub = 'COMPLETE'

        # We're done. Return the result.
        return result
