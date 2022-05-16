from ..stpipe import Step
from . import ami_average

__all__ = ["AmiAverageStep"]


class AmiAverageStep(Step):
    """
    AmiAverageStep: Averages LG results for multiple NIRISS AMI mode exposures
    """

    class_alias = "ami_average"

    spec = """
    """

    def flatten_input(self, input_items):
        """Remove any nested list/tuple structure and return generator
        to provide iterable simple list with no nested structure.
        """
        for item in input_items:
            if isinstance(item, (list, tuple)):
                yield from self.flatten_input(item)
            else:
                yield item

    def process(self, *input_list):
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

        # Input may be a simple list if run in interactive environment,
        # but processing from command line wraps list of inputs in a tuple.
        # Flatten this object into a simple list.
        input_list = list(self.flatten_input(input_list))

        # Call the LG average routine for the list of inputs
        result = ami_average.average_LG(input_list)

        result.meta.cal_step.ami_average = 'COMPLETE'

        return result
