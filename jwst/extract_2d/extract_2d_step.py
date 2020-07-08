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
        extract_orders = int_list(default=None)  # list of orders to extract
        extract_height =  integer(default=None)  # extraction height in pixels
        grism_objects = list(default=None)  # list of grism objects to use
        mmag_extract = float(default=99.)  # minimum abmag to extract
    """

    reference_file_types = ['wavelengthrange']

    def process(self, input_model, *args, **kwargs):
        reffile = self.get_reference_file(input_model, 'wavelengthrange')
        reference_file_name = reffile if reffile else ""

        with datamodels.open(input_model) as dm:
            output_model = extract_2d.extract2d(dm, self.slit_name,
                                                reference_file=reference_file_name,
                                                extract_orders=self.extract_orders,
                                                grism_objects=self.grism_objects,
                                                extract_height=self.extract_height)

        return output_model
