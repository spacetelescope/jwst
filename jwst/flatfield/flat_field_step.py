#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import flat_field

# For the following types of data, it is OK -- and in some cases
# required -- for the extract_2d step to have been run.  For all
# other types of data, the extract_2d step must not have been run.
EXTRACT_2D_IS_OK = ["NRS_LAMP", "NRS_BRIGHTOBJ", "NRS_FIXEDSLIT",
                    "NRS_MSASPEC"]

# NIRSpec imaging types (see exp_type2transform in assign_wcs/nirspec.py).
IMAGING_MODES = ["NRS_IMAGE", "NRS_FOCUS",
                 "NRS_TACQ", "NRS_BOTA", "NRS_TASLIT",
                 "NRS_CONFIRM", "NRS_TACONFIRM",
                 "NRS_MIMF"]

class FlatFieldStep(Step):
    """
    FlatFieldStep: Flat-field a science image using a flatfield reference image.
    """

    spec = """
        # Suffix for optional output file for interpolated flat fields.
        # Note that this is only used for NIRSpec spectrographic data.
        flat_suffix = string(default=None)
    """

    reference_file_types = ["flat", "fflat", "sflat", "dflat"]

    def process(self, input):

        if self.flat_suffix is not None:
            if self.flat_suffix == "None" or len(self.flat_suffix) == 0:
                self.flat_suffix = None

        input_model = datamodels.open(input)
        exposure_type = input_model.meta.exposure.type.upper()

        # Figure out what kind of input data model is in use.
        if isinstance(input_model, datamodels.CubeModel):
            # multi-integration dataset
            self.log.debug('Input is a CubeModel')
        elif isinstance(input_model, datamodels.ImageModel):
            self.log.debug('Input is an ImageModel')
        elif isinstance(input_model, datamodels.MultiSlitModel):
            self.log.debug('Input is a MultiSlitModel')

        # Check whether extract_2d has been run.
        if (input_model.meta.cal_step.extract_2d == 'COMPLETE' and
            not exposure_type in EXTRACT_2D_IS_OK):
            self.log.warning("The extract_2d step has been run, but for "
                             "%s data it should not have been run, so ...",
                             exposure_type)
            self.log.warning("flat fielding will be skipped.")
            result = input_model.copy()
            result.meta.cal_step.flat_field = "SKIPPED"
            input_model.close()
            return result

        # NIRSpec spectrographic mode?
        if input_model.meta.instrument.name.upper() == "NIRSPEC":
            if exposure_type in IMAGING_MODES:
                is_NRS_spectrographic = False
            else:
                is_NRS_spectrographic = True
        else:
            is_NRS_spectrographic = False

        # Retrieve the reference file name or names
        if is_NRS_spectrographic:
            self.flat_filename = 'N/A'
            self.f_flat_filename = self.get_reference_file(input_model,
                                        'fflat')
            self.s_flat_filename = self.get_reference_file(input_model,
                                        'sflat')
            self.d_flat_filename = self.get_reference_file(input_model,
                                        'dflat')
            self.log.debug('Using FFLAT reference file: %s',
                           self.f_flat_filename)
            self.log.debug('Using SFLAT reference file: %s',
                           self.s_flat_filename)
            self.log.debug('Using DFLAT reference file: %s',
                           self.d_flat_filename)
        else:
            self.flat_filename = self.get_reference_file(input_model, 'flat')
            self.f_flat_filename = 'N/A'
            self.s_flat_filename = 'N/A'
            self.d_flat_filename = 'N/A'
            self.log.debug('Using FLAT reference file: %s', self.flat_filename)

        # Check for a valid reference file
        missing = False
        if is_NRS_spectrographic:
            if (self.f_flat_filename == 'N/A' or
                self.s_flat_filename == 'N/A' or
                self.d_flat_filename == 'N/A'):
                self.log.warning('One or more flat-field reference files'
                                 ' were missing')
                missing = True
        else:
            if self.flat_filename == 'N/A':
                self.log.warning('No FLAT reference file found')
                missing = True
        if missing:
            self.log.warning('Flat-field step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.flat_field = 'SKIPPED'
            input_model.close()
            return result

        # Find out what model to use for the flat field reference file(s).
        if is_NRS_spectrographic:
            flat_model = None
            if exposure_type == "NRS_MSASPEC":
                f_flat_model = \
                    datamodels.NirspecQuadFlatModel(self.f_flat_filename)
            else:
                f_flat_model = \
                    datamodels.NirspecFlatModel(self.f_flat_filename)
            s_flat_model = datamodels.NirspecFlatModel(self.s_flat_filename)
            d_flat_model = datamodels.NirspecFlatModel(self.d_flat_filename)
        else:
            self.log.debug('Opening flat as FlatModel')
            flat_model = datamodels.FlatModel(self.flat_filename)
            f_flat_model = None
            s_flat_model = None
            d_flat_model = None
            self.flat_suffix = None

        # Do the flat-field correction
        (output_model, interpolated_flats) = \
                flat_field.do_correction(input_model, flat_model,
                                         f_flat_model, s_flat_model,
                                         d_flat_model, self.flat_suffix)

        # Close the inputs
        input_model.close()
        if is_NRS_spectrographic:
            f_flat_model.close()
            s_flat_model.close()
            d_flat_model.close()
        else:
            flat_model.close()

        if interpolated_flats is not None:
            self.log.info("Writing interpolated flat fields.")
            self.save_model(interpolated_flats, self.flat_suffix)
            interpolated_flats.close()

        return output_model

if __name__ == '__main__':
    cmdline.step_script(flat_field_step)
