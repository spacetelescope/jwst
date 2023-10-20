#! /usr/bin/env python

from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
from ..stpipe import Step
import numpy as np

from . import spectral_leak

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

        ch1b = None
        ch3a = None
        ich3a = None
        ch1b_wave = 6.0
        ch3a_wave = 12.0

        with datamodels.open(input) as input_model:
            if isinstance(input_model, ModelContainer):
                self.log.info('Input is a ModelContainer')
                for i, x1d in enumerate(input_model):
                    channel = x1d.meta.instrument.channel
                    band = x1d.meta.instrument.band
                    srctype = x1d.spec[0].source_type
                    if srctype == 'EXTENDED':
                        self.log.info('No spectral leak correction for extended source data')
                        return input
                    # search x1d containing CH 1 B
                    if  '1' in channel and 'MEDIUM' in band:
                        print('found ch 1B')
                        ch1b = x1d 
                    elif '1' in channel and 'MULTIPLE' in band:
                        # read in the wavelength array and see
                        # if it covers ch1b_wave
                        wave = x1d.spec[0].spec_table.WAVELENGTH
                        if np.min(wave) < ch1b_wave and np.max(wave) > ch1b_wave:
                            print('found ch 1B from wavelength')
                            ch1b = x1d
                     # search x1d containing CH 3 A
                    if  '3' in channel and 'SHORT' in band:
                        print('found ch 3A')
                        ch3a = x1d
                        ich3a = i # store the datamodel to update later 
                    elif '#' in channel and 'MULTIPLE' in band:
                        # read in the wavelength array and see
                        # if it covers ch3a_wave
                        wave = x1d.spec[0].spec_table.WAVELENGTH
                        if np.min(wave) < ch3a_wave and np.max(wave) > ch3a_wave:
                            print('found ch 3A from wavelength')
                            ch3a = x1d
                            ich3a = i # store the datamodel to update later
        if ch1b is not None and ch3a is not None:
            sp_leak_ref = self.get_reference_file(input[0], 'mrsptcorr')
            corrected_3a = spectral_leak.do_correction(sp_leak_ref, ch1b, ch3a)
            input[ich3a].spec[0].spec_table.FLUX = corrected_3a
            input[ich3a].meta.cal_step.spectral_leak = 'COMPLETE'
            input[ich3a].save('test.fits')
            print(input.meta.ref_file.name)
            print(input[ich3a].meta.instrument.channel)
            print(input[ich3a].meta.instrument.band)

        return input

