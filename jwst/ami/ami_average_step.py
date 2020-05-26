from ..stpipe import Step

from . import ami_average

__all__ = ["AmiAverageStep"]


class AmiAverageStep(Step):
    """
    AmiAverageStep: Averages LG results for multiple NIRISS AMI mode exposures
    """

    spec = """
    """

    def process(self, input_list):
        """
        Averages the results of LG analysis for a set of multiple NIRISS AMI
        mode exposures.

        Parameters
        ----------
        input_list: list
            input file names

        Returns
        -------
        result: AmiLgModel object
            Averaged AMI data model
        """

        # Call the LG average routine for the list of inputs
        result = ami_average.average_LG(input_list)

        result.meta.cal_step.ami_average = 'COMPLETE'

        return result
