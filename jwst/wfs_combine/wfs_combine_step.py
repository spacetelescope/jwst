#! /usr/bin/env python

from ..stpipe import Step, cmdline
from . import wfs_combine
import json
import os

class WfsCombineStep(Step):

    """
    This step combines pairs of dithered PSF images
    """

    spec = """
        do_refine = boolean(default=False)
    """

    def process(self, input_table):

        asn_table = json.load(open(input_table, 'r'))
        num_sets = len(asn_table['products'])

        self.log.info('Using input table: %s', input_table)
        self.log.info('The number of pairs of input files: %g', num_sets)

        # Process each pair of input images listed in the association table
        for which_set in range(num_sets):
            infile_1 = asn_table['products'][which_set]['members'][0]['expname']
            infile_2 = asn_table['products'][which_set]['members'][1]['expname']
            outfile = asn_table['products'][which_set]['name']

            # Construct the full output file name
            outfile = self.make_output_path(
                None, basepath=outfile, suffix='wfscmb'
            )

            wfs = wfs_combine.DataSet(
                infile_1, infile_2, outfile, self.do_refine
            )

            output_model = wfs.do_all()
            output_model.meta.cal_step.wfs_combine = 'COMPLETE'
            self.save_model(
                output_model, 'wfscmb', output_file=outfile
            )

        return None


def mk_prodname(output_dir, filename, suffix):

    if output_dir is not None:
        dirname, filename = os.path.split(filename)
        filename = os.path.join(output_dir, filename)

    base, ext = os.path.splitext(filename)
    if len(ext) == 0:
        ext = ".fits"
    return base + '_' + suffix + ext
