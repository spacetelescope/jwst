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
        transmission_level = integer(default=50) # transmission level to use to define slice locations
    """
    
    reference_file_types = ['fringefreq','regions']

    def process(self, input_data):

        self.suffix = 'residual_fringe'
        
        with datamodels.open(input_data) as input_models:

            self.input_models = input_models
            if not isinstance(self.input_models, datamodels.ModelContainer):
                self.input_container = False
            else:
                self.input_container = True

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

            if self.input_container:
                exptype = self.input_models[0].meta.exposure.type
            else:
                exptype = self.input_models.meta.exposure.type

            print('exptype',exptype)
            if exptype != 'MIR_MRS':
                self.log(" Residual Fringe correction is only for MIRI MRS data")
                self.log.error( "Unsupported ", f"exposure type: {exptype}")

                if self.input_container:
                    for model in self.input_models:
                        model.meta.cal_step.outlier_detection = "SKIPPED"
                else:
                    self.input_models.meta.cal_step.outlier_detection = "SKIPPED"
                self.skip = True

                return self.input_models

        # loop over each model
        # read in reference file for model
        # form ResidualFringeCorrection
        # correct
        for model in self.input_models:
            
            # Open the residual fringe reference file
            self.residual_fringe_filename = self.get_reference_file(model,
                                                                    'fringefreq')
            self.log.info('Using Residual FRINGE reference file: %s',
                          self.residual_fringe_filename)

            # Open the regions reference file
            self.regions_filename = self.get_reference_file(model,'regions')
            self.log.info('Using MRS regrions reference file: %s',
                          self.regions_filename)

            # Check for a valid reference files. If they are not found skip step 
            if self.residual_fringe_filename == 'N/A' or self.regions_filename =='N/A':
                
                if self.residual_fringe_filename == 'N/A':
                    self.log.warning('No Residual FRINGE reference file found')
                    self.log.warning('Residual Fringe step will be skipped')

                if self.regions_filename == 'N/A':
                    self.log.warning('No MRS regions reference file found')
                    self.log.warning('Residual Fringe step will be skipped')

                if self.input_container:
                    for model in self.input_models:
                        model.meta.cal_step.outlier_detection = "SKIPPED"
                else:
                    self.input_models.meta.cal_step.residual_fringe = "SKIPPED"
                self.skip = True
                return self.input_models

            # Do the correction
            rfc= residual_fringe.ResidualFringeCorrection(model,
                                                          self.residual_fringe_filename,
                                                          self.regions_filename,
                                                          self.transmission_level)
            rfc.do_correction()






        return output_model


    def make_output_path(ignored, idx=None):
                output_path = self.make_output_path(
                    basepath=base_filename, suffix='rf', idx=idx,
                    component_format='_{asn_id}_{idx}'
                )
                return output_path
