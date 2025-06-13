#! /usr/bin/env python

from stdatamodels.jwst import datamodels
from jwst.datamodels import ModelContainer
from jwst.stpipe import Step
import numpy as np

from . import spectral_leak

__all__ = ["SpectralLeakStep"]


class SpectralLeakStep(Step):
    """
    Apply a spectral leak correction to the Channel 3A of MIRI MRS data.

    The MIRI MRS has a spectral leak in which 6 micron light leaks into the
    12 micron channel.  This step applies a correction to the 12 micron channel.
    """

    class_alias = "spectral_leak"

    reference_file_types = ["mrsptcorr"]

    def process(self, input_data):
        """
        Apply a spectral leak  correction to MIRI MRS Channel 3A data.

        Parameters
        ----------
        input_data : ~jwst.datamodels.ModelContainer
            Container of models containing 1-D extracted spectra

        Returns
        -------
        JWST DataModel
            The corrected data model. This will be the input model if the step is skipped,
            otherwise it will be a corrected 1D extracted spectrum that contains
            the MRS channel 3B range.
        """
        ch1b = None
        ch3a = None
        ich3a = None
        ch1b_wave = 6.0
        ch3a_wave = 12.0

        with datamodels.open(input_data) as input_model:
            if isinstance(input_model, ModelContainer):
                result = input_model.copy()  # copy input to return
                # Retrieve the reference parameters for this type of data
                sp_leak_ref = self.get_reference_file(input_model[0], "mrsptcorr")

                for i, x1d in enumerate(input_model):  # input_model is a Model Container
                    # check that we have the correct type of data
                    if isinstance(x1d, datamodels.MultiSpecModel):
                        self.log.debug(" Data is MIRI MRS MultiSpecModel data")
                    elif isinstance(x1d, datamodels.MRSMultiSpecModel):
                        self.log.debug(" Data is  MIRI MRS MRSMultiSpecModel data")
                    else:
                        self.log.warning(
                            "Data sent to spectral_leak step is not an extracted spectrum. "
                            f" It is  {type(x1d)}."
                        )
                        for r in result:
                            r.meta.cal_step.spectral_leak = "SKIPPED"
                        return result
                    channel = x1d.meta.instrument.channel
                    band = x1d.meta.instrument.band
                    srctype = x1d.spec[0].source_type
                    if srctype == "EXTENDED":
                        self.log.warning("No spectral leak correction for extended source data")
                        for r in result:
                            r.meta.cal_step.spectral_leak = "SKIPPED"
                        return result
                    # search x1d containing CH 1 B
                    if "1" in channel and "MEDIUM" in band:
                        self.log.info("Found CH 1B in input data")
                        ch1b = x1d
                    elif "1" in channel and "MULTIPLE" in band:
                        # read in the wavelength array and see
                        # if it covers ch1b_wave
                        wave = x1d.spec[0].spec_table.WAVELENGTH
                        if np.min(wave) < ch1b_wave and np.max(wave) > ch1b_wave:
                            self.log.info("Found CH 1B in the input data")
                            ch1b = x1d
                    # search x1d containing CH 3 A
                    if "3" in channel and "SHORT" in band:
                        self.log.info("Found CH 3A in the input data")
                        ch3a = x1d
                        ich3a = i  # store the datamodel # to update later
                    elif "3" in channel and "MULTIPLE" in band:
                        # read in the wavelength array and see
                        # if it covers ch3a_wave
                        wave = x1d.spec[0].spec_table.WAVELENGTH
                        if np.min(wave) < ch3a_wave and np.max(wave) > ch3a_wave:
                            self.log.info("Found CH 3A in the input data")
                            ch3a = x1d
                            ich3a = i  # store the datamodel to update later

                # done looping over data now if 1B and 3A data exists make a correction
                # update result and return
                if ch1b is not None and ch3a is not None:
                    corrected_3a, corrected_3a_rf = spectral_leak.do_correction(
                        sp_leak_ref, ch1b, ch3a
                    )
                    result[ich3a].spec[0].spec_table.FLUX = corrected_3a
                    if corrected_3a_rf is not None:
                        result[ich3a].spec[0].spec_table.RF_FLUX = corrected_3a_rf
                    result[ich3a].meta.cal_step.spectral_leak = "COMPLETE"
                    return result
                else:
                    for r in result:
                        r.meta.cal_step.spectral_leak = "SKIPPED"
                    self.log.warning("CH1B and CH3A were not found. No spectral leak correction")
                    return result

            else:
                self.log.warning(
                    "Data sent to spectral_leak step is not a ModelContainer."
                    "It is {type(input_model)}."
                )
                self.log.warning("Step is skipped")
                return input_data
