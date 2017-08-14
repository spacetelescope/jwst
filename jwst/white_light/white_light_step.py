#! /usr/bin/env python

from ..stpipe import Step, cmdline
from .. import datamodels
from ..source_catalog.source_catalog import _replace_suffix_ext
from .white_light import white_light


class WhiteLightStep(Step):
    """
    WhiteLightStep: Computes integrated flux as a function of time for a
    multi-integration spectroscopic observation.
    """

    spec = """
    """

    def process(self, input):

        # Load the input
        with datamodels.open(input) as input_model:

            # Call the white light curve generation routine
            result = white_light(input_model)

            # Write the output catalog
            old_suffixes = ['x1dints']
            output_dir = self.search_attr('output_dir')
            cat_filepath = _replace_suffix_ext(
                input_model.meta.filename, old_suffixes, 'whtlt',
                output_ext='ecsv', output_dir=output_dir)
            result.write(cat_filepath, format='ascii.ecsv', overwrite=True)

        return
