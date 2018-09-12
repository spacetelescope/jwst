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
        slit_name = string(default=None)
        apply_wavecorr = boolean(default=True)
        extract_orders = int_list(default=None)  # list of orders to extract
        extract_height =  integer(default=None)  # extraction height in pixels
        grism_objects = list(default=None)  # list of grism objects to use
    """

    reference_file_types = ['wavecorr', 'wavelengthrange']

    def process(self, input_model, *args, **kwargs):
        reference_file_names = {}
        for reftype in self.reference_file_types:
            reffile = self.get_reference_file(input_model, reftype)
            reference_file_names[reftype] = reffile if reffile else ""

        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(dm, self.slit_name, self.apply_wavecorr,
                                                reference_files=reference_file_names,
                                                extract_orders=self.extract_orders,
                                                grism_objects=self.grism_objects,
                                                extract_height=self.extract_height)

        return output_model
