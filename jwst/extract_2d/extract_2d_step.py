#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from . import extract_2d


class Extract2dStep(Step):
    """
    This Step performs a 2D extraction of spectra.
    """

    spec = """
        which_subarray = string(default = None)
        apply_wavecorr = boolean(default=True)
    """

    reference_file_types = ['wavecorr']


    def process(self, input_model, *args, **kwargs):
        reffile = self.get_reference_file(input_model, "wavecorr")
        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(dm, self.which_subarray,
                                                self.apply_wavecorr, reffile=reffile)

        return output_model


if __name__ == '__main__':
    cmdline.step_script(Extract2dStep)
