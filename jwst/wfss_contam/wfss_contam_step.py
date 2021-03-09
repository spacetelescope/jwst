#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import wfss_contam


__all__ = ["WfssContamStep"]


class WfssContamStep(Step):
    """
    This Step performs contamination correction of WFSS spectra.
    """

    spec = """
    """

    reference_file_types = ['wavelengthrange']


    def process(self, input_model, *args, **kwargs):

        with datamodels.open(input_model) as dm:

            # Get the wavelengthrange ref file
            waverange_ref = self.get_reference_file(dm, 'wavelengthrange')
            self.log.info(f'Using WAVELENGTHRANGE reference file {waverange_ref}')
            waverange_model = datamodels.WavelengthrangeModel(waverange_ref)

            output_model = wfss_contam.contam_corr(dm, waverange_model)

        return output_model
