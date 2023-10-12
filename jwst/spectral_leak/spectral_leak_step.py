#! /usr/bin/env python

from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
from ..stpipe import Step


#from . import spectral_leak

__all__ = ["SpectralLeakStep"]


class SpectralLeakStep(Step):
    """
    The MIRI MRS has a spectral leak in which 6 micron light leaks into the 
    12 micron channel.  This step applies a correction to the 12 micron channel.
    """

    class_alias = "spectral_leak"

    reference_file_types = ['mrsptcorr']

    def process(self, input):
        """Execute the step.

        Parameters
        ----------
        input : container of models containing 1-D extracted spectra

        Returns
        -------
        JWST DataModel
            This will be `input` if the step was skipped; otherwise,
            it will be a corrected 1-D extracted spectra that contains
            the  3B MRS range. 
        """

        
        # Open the input and figure out what type of model it is
        input_model = datamodels.open(input)

        found_1b = None
        found_3a = None
        with datamodels.open(input) as input_model:
            if isinstance(input_model, ModelContainer):
                self.log.debug('Input is a ModelContainer')
                channel = input_model.meta.instrument.channel
                band = input_model.meta.instrument.band
                srctype = input_model.meta.target.source_type
                print(channel, band, srctype)
                
                #apcorr_ref = (
                #        self.get_reference_file(input_model[0], 'apcorr')
                    
            #return r

