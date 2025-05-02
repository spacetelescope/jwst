from stdatamodels.jwst import datamodels

from ..stpipe import Step

from . import stack_refs

__all__ = ["StackRefsStep"]


class StackRefsStep(Step):
    """
    Stack multiple PSF reference exposures into a single CubeModel.

    Result used by subsequent coronagraphic steps.
    """

    class_alias = "stack_refs"

    spec = """
    """  # noqa: E501

    def process(self, input_files):
        """
        Execute the StackRefsStep calibration step.

        Parameters
        ----------
        input_files : str
            Input science exposures

        Returns
        -------
        output_model : DataModel
            PSF reference exposures stacked into a CubeModel
        """
        # Open the inputs
        with datamodels.open(input_files) as input_models:
            # Call the stacking routine
            output_model = stack_refs.make_cube(input_models)

            output_model.meta.cal_step.stack_psfs = "COMPLETE"

        return output_model
