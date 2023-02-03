#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import residual_fringe
from functools import partial

__all__ = ["ResidualFringeStep"]


class ResidualFringeStep(Step):
    """
    ResidualFringeStep: Apply residual fringe correction to a science image
    using parameters in the residual fringe reference file.

    Parameters
    ----------
    input_data : asn file or single file
    """

    class_alias = 'residual_fringe'

    spec = """
        skip = boolean(default=True)
        save_intermediate_results  = boolean(default = False)
        search_output_file = boolean(default = False)
        ignore_region_min = list(default = None)
        ignore_region_max = list(default = None)
        suffix = string(default = 'residual_fringe')
    """

    reference_file_types = ['fringefreq', 'regions']

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

        if isinstance(input, datamodels.IFUImageModel):
            exptype = input.meta.exposure.type
        else:
            raise TypeError("Failed to process file type {}".format(type(input)))

        # Setup output path naming if associations are involved.
        asn_id = None
        try:
            asn_id = self.input.meta.asn_table.asn_id
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
            input.meta.cal_step.residual_fringe = "SKIPPED"
            return input

        # 1. set up the reference files
        # 2. correct the  model
        # 3. return from step

        self.residual_fringe_filename = self.get_reference_file(input, 'fringefreq')
        self.log.info('Using FRINGEFREQ reference file:{}'.
                      format(self.residual_fringe_filename))

        # set up regions reference file
        self.regions_filename = self.get_reference_file(input, 'regions')
        self.log.info('Using MRS regions reference file: {}'.
                      format(self.regions_filename))

        # Check for a valid reference files. If they are not found skip step
        if self.residual_fringe_filename == 'N/A' or self.regions_filename == 'N/A':
            if self.residual_fringe_filename == 'N/A':
                self.log.warning('No FRINGEFREQ reference file found')
                self.log.warning('Residual Fringe step will be skipped')

            if self.regions_filename == 'N/A':
                self.log.warning('No MRS regions reference file found')
                self.log.warning('Residual Fringe step will be skipped')

            input.meta.cal_step.residual_fringe = "SKIPPED"
            return input

        # Do the correction
        rfc = residual_fringe.ResidualFringeCorrection(input,
                                                       self.residual_fringe_filename,
                                                       self.regions_filename,
                                                       self.ignore_regions,
                                                       **pars)
        result = rfc.do_correction()
        result.meta.cal_step.residual_fringe = 'COMPLETE'
        return result
