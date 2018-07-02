#! /usr/bin/env python

from ..stpipe import Step
from .. import datamodels
from . import extract_2d


__all__ = ["Extract2dStep"]


class Extract2dStep(Step):
    """
    This Step performs a 2D extraction of spectra.
    """

    spec = """
        slit_name = string(default = None)
        apply_wavecorr = boolean(default=True)
        grism_objects = list(default=list())
    """

    reference_file_types = ['wavecorr', 'wavelengthrange']

    def process(self, input_model, *args, **kwargs):
        reference_file_names = {}
        for reftype in self.reference_file_types:
            reffile = self.get_reference_file(input_model, reftype)
            reference_file_names[reftype] = reffile if reffile else ""

        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(dm, self.slit_name, self.apply_wavecorr,
                                                reference_files=reference_file_names)

        return output_model
