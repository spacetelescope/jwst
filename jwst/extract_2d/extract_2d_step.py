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

    reference_file_types = ['wavecorr', 'wavelengthrange']

    def process(self, input_model, *args, **kwargs):
        reference_file_names = {}
        for reftype in self.reference_file_types:
            reffile = self.get_reference_file(input_model, reftype)
            reference_file_names[reftype] = reffile if reffile else ""

        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(dm, self.which_subarray, self.apply_wavecorr,
                                                reference_files=reference_file_names)

        return output_model

