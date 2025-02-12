#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from ..stpipe import Step
from . import straylight

__all__ = ["StraylightStep"]


class StraylightStep (Step):
    """
    StraylightStep: Performs straylight correction image using a Mask file.
    """

    class_alias = "straylight"

    spec = """
        clean_showers = boolean(default=False)  # Clean up straylight from residual cosmic ray showers
        shower_plane = integer(default=3)  # Throughput plane for identifying inter-slice regions
        shower_x_stddev = float(default=18) # X standard deviation for shower model
        shower_y_stddev = float(default=5) # Y standard deviation for shower model
        shower_low_reject = float(default=0.1) # Low percentile of pixels to reject
        shower_high_reject = float(default=99.9) # High percentile of pixels to reject
    """  # noqa: E501

    reference_file_types = ['mrsxartcorr', 'regions']

    def process(self, input):

        with datamodels.open(input) as input_model:
            # Set up the output result
            result = input_model.copy()

            # check the data is an IFUImageModel (not TSO)
            if isinstance(input_model, (datamodels.ImageModel, datamodels.IFUImageModel)):
                # Check for a valid mrsxartcorr reference file
                self.straylight_name = self.get_reference_file(input_model, 'mrsxartcorr')

                if self.straylight_name == 'N/A':
                    self.log.warning('No MRSXARTCORR reference file found')
                    self.log.warning('Straylight step will be skipped')
                    result.meta.cal_step.straylight = 'SKIPPED'
                    return result

                self.log.info('Using mrsxartcorr reference file %s', self.straylight_name)

                modelpars = datamodels.MirMrsXArtCorrModel(self.straylight_name)

                # Apply the cross artifact correction
                result = straylight.correct_xartifact(result, modelpars)

                modelpars.close()

                # Apply the cosmic ray droplets correction if desired
                if self.clean_showers:
                    self.regions_name = self.get_reference_file(input_model, 'regions')
                    with datamodels.RegionsModel(self.regions_name) as f:
                        allregions = f.regions
                        result = straylight.clean_showers(result, allregions, self.shower_plane,
                                                          self.shower_x_stddev, self.shower_y_stddev,
                                                          self.shower_low_reject, self.shower_high_reject)

                result.meta.cal_step.straylight = 'COMPLETE'

            else:
                if isinstance(input_model, (datamodels.ImageModel, datamodels.IFUImageModel)) is False:
                    self.log.warning('Straylight correction not defined for datatype %s',
                                     input_model)
                self.log.warning('Straylight step will be skipped')
                result.meta.cal_step.straylight = 'SKIPPED'

        return result


class ErrorNoAssignWCS(Exception):
    pass
