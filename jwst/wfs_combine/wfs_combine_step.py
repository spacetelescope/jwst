#! /usr/bin/env python

from ..stpipe import Step
from . import wfs_combine

__all__ = ["WfsCombineStep"]


class WfsCombineStep(Step):

    """
    This step combines pairs of dithered PSF images
    """

    spec = """
        do_refine = boolean(default=False)
    """

    def process(self, input_table):

        asn_table = self.load_as_level3_asn(input_table)
        num_sets = len(asn_table['products'])

        self.log.info('Using input table: %s', input_table)
        self.log.info('The number of pairs of input files: %g', num_sets)

        # Process each pair of input images listed in the association table
        for which_set in range(num_sets):
            infile_1 = asn_table['products'][which_set]['members'][0]['expname']
            infile_2 = asn_table['products'][which_set]['members'][1]['expname']
            outfile = asn_table['products'][which_set]['name']

            wfs = wfs_combine.DataSet(
                infile_1, infile_2, outfile, self.do_refine
            )

            output_model = wfs.do_all()
            output_model.meta.cal_step.wfs_combine = 'COMPLETE'
            self.save_model(
                output_model, suffix='wfscmb', output_file=outfile, format=False
            )

        return None
