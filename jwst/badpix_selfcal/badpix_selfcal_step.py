

from stpipe import Step
from . import badpix_selfcal


class BadpixSelfcalStep(Step):
    """
    BadpixSelfcalStep: 
    """

    class_alias = "badpix_selfcal"

    # Define a suffix for optional saved output
    bkg_suffix = 'badpix_selfcal'

    spec = """
    flagpercentile = float(default=99.9)  # percentile level above which to flag
    skip = boolean(default=True)
    """

    def process(self, input):
        """
        Flag CRs in the DQ array of a JWST exposure

        Parameters
        ----------
        input: JWST data model or association
            input science data to be corrected

        Returns
        -------
        output: JWST data model or association
            data model with CRs flagged
        """
        # Do the work here
        return result