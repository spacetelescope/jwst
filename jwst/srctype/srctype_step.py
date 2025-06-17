#! /usr/bin/env python
from stdatamodels.jwst import datamodels

from jwst.stpipe import Step
from .srctype import set_source_type

__all__ = ["SourceTypeStep"]


class SourceTypeStep(Step):
    """
    Select and set a source type based on various inputs.

    The source type is used in later calibrations to determine the appropriate
    methods to use. Input comes from either the SRCTYAPT keyword value, which
    is populated from user info in the APT, or the NIRSpec MSA planning tool.
    The source type can be also specified on the command line for exposures
    containing a single pre-defined target.
    """

    class_alias = "srctype"

    spec = """
        source_type = option('POINT','EXTENDED', default=None)  # user-supplied source type
    """

    def process(self, step_input):
        """
        Determine the source type.

        Parameters
        ----------
        step_input : str, IFUImageModel, MultiSlitModel, or SlitModel
            Either the path to the file or the science data model to be corrected.

        Returns
        -------
        datamodel : IFUImageModel, MultiSlitModel, or SlitModel
            Data model with keyword “SRCTYPE” populated with either “POINT” or “EXTENDED”.
        """
        if self.source_type is not None:
            self.source_type = self.source_type.upper()

        source_type = self.source_type  # retrieve command line override

        input_model = datamodels.open(step_input)

        # Call the source selection routine
        result = set_source_type(input_model, source_type)

        # Set the step status in the output model
        if result is None:
            result = input_model
            result.meta.cal_step.srctype = "SKIPPED"
        else:
            result.meta.cal_step.srctype = "COMPLETE"

        return result
