#! /usr/bin/env python
from ..stpipe import Step
from .. import datamodels
from . import residual_fringe
from functools import partial

__all__ = ["ResidualFringeStep"]


class ResidualFringeStep(Step):
    """
    ResidualFringeStep: Apply residual fringe correction to a science image using a residual fringe
    reference image.

    Parameters
    ----------
    input_data : asn file or ModelContainer
        Single filename association table or a datamodels.ModelContainer
    """

    spec = """
        save_intermediate_results  = boolean(default = False)
        search_output_file = boolean(default = False)
        ignore_region_min = list(default = None)
        ignore_region_max = list(default = None)
        suffix = string(default = 'residual_fringe')
    """

    reference_file_types = ['fringefreq', 'regions']

    class_alias = 'residual_fringe'

    def process(self, input):

        self.transmission_level = 80  # sets the transmission level to use in the regions file
        # 80% is what other steps use.

        # set up the dictionary to ignore wavelength regions in the residual fringe correction
        ignore_regions = {}
        ignore_regions['num'] = 0
        ignore_regions['min'] = []
        ignore_regions['max'] = []
        if self.ignore_region_min is not None:
            for region in self.ignore_region_min:
                ignore_regions['min'].append(float(region))

        min_num = len(ignore_regions['min'])

        if self.ignore_region_max is not None:
            for region in self.ignore_region_max:
                ignore_regions['max'].append(float(region))
        max_num = len(ignore_regions['max'])

        if max_num != min_num:
            self.log.error("Number of minimum and maximum wavelengths to ignore are not the same")
            raise ValueError("Number of ignore_region_min does not match ignore_region_max")

        ignore_regions['num'] = min_num

        if min_num > 0:
            self.log.info('Ignoring {} wavelength regions'.format(min_num))

        self.ignore_regions = ignore_regions

        input = datamodels.open(input)

        # If single file, wrap in a ModelContainer
        if isinstance(input, datamodels.IFUImageModel):
            input_models = datamodels.ModelContainer([input])
            self.input_container = False
            exptype = input.meta.exposure.type
        elif isinstance(input, datamodels.ModelContainer):
            input_models = input
            self.input_container = True
            exptype = input[0].meta.exposure.type
        else:
            raise TypeError("Failed to process file type {}".format(type(input)))

        # Setup output path naming if associations are involved.
        asn_id = None
        try:
            asn_id = self.input_models.meta.asn_table.asn_id
        except (AttributeError, KeyError):
            pass
        if asn_id is None:
            asn_id = self.search_attr('asn_id')
        if asn_id is not None:
            _make_output_path = self.search_attr(
                '_make_output_path', parent_first=True
            )

            self._make_output_path = partial(
                _make_output_path,
                asn_id=asn_id
            )

        # Set up residual fringe correction parameters
        pars = {
            'transmission_level': self.transmission_level,
            'save_intermediate_results': self.save_intermediate_results,
            'make_output_path': self.make_output_path
        }

        if exptype != 'MIR_MRS':
            self.log(" Residual Fringe correction is only for MIRI MRS data")
            self.log.error("Unsupported ", f"exposure type: {exptype}")

            if self.input_container:
                for model in input_models:
                    model.meta.cal_step.outlier_detection = "SKIPPED"
                else:
                    input_models.meta.cal_step.outlier_detection = "SKIPPED"
                self.skip = True
                return input_models

        # loop over each model
        # 1. read in reference file for model
        # 2. correct each model
        # 3. append corrected data to output_models - to return from step
        self.output_models = datamodels.ModelContainer()
        for model in input_models:
            # Open the residual fringe reference file
            self.residual_fringe_filename = self.get_reference_file(model,
                                                                    'fringefreq')
            self.log.info('Using Residual FRINGE reference file:{}'.
                          format(self.residual_fringe_filename))

            # Open the regions reference file
            self.regions_filename = self.get_reference_file(model, 'regions')
            self.log.info('Using MRS regions reference file: {}'.
                          format(self.regions_filename))

            # Check for a valid reference files. If they are not found skip step
            if self.residual_fringe_filename == 'N/A' or self.regions_filename == 'N/A':

                if self.residual_fringe_filename == 'N/A':
                    self.log.warning('No Residual FRINGE reference file found')
                    self.log.warning('Residual Fringe step will be skipped')

                if self.regions_filename == 'N/A':
                    self.log.warning('No MRS regions reference file found')
                    self.log.warning('Residual Fringe step will be skipped')

                if self.input_container:
                    for model in self.input_models:
                        model.meta.cal_step.outlier_detection = "SKIPPED"
                        self.output_models.append(model)
                else:
                    input.meta.cal_step.residual_fringe = "SKIPPED"
                    self.output_models.append(input)
                    self.skip = True
                return self.output_models

            # Do the correction
            rfc = residual_fringe.ResidualFringeCorrection(model,
                                                           self.residual_fringe_filename,
                                                           self.regions_filename,
                                                           self.ignore_regions,
                                                           **pars)
            result = rfc.do_correction()
            result.meta.cal_step.residual_fringe = 'COMPLETE'
            self.output_models.append(result)

        return self.output_models
