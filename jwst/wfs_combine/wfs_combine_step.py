import os

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

        # Load the input ASN table
        asn_table = self.load_as_level3_asn(input_table)
        num_sets = len(asn_table['products'])

        self.log.info('Using input table: %s', input_table)
        self.log.info('The number of pairs of input files: %g', num_sets)

        # Process each pair of input images listed in the association table
        for which_set in asn_table['products']:

            # Get the list of science members in this pair
            science_members = [
                member
                for member in which_set['members']
                if member['exptype'].lower() == 'science'
            ]
            infile_1 = science_members[0]['expname']
            infile_2 = science_members[1]['expname']
            outfile = which_set['name']

            # Create the step instance
            wfs = wfs_combine.DataSet(
                infile_1, infile_2, outfile, self.do_refine
            )

            # Do the processing
            output_model = wfs.do_all()

            # The DataSet class does not close its resources.  Do that here.
            wfs.input_1.close()
            wfs.input_2.close()

            # Update necessary meta info in the output
            output_model.meta.cal_step.wfs_combine = 'COMPLETE'
            output_model.meta.asn.pool_name = asn_table['asn_pool']
            output_model.meta.asn.table_name = os.path.basename(input_table)

            # Save the output file
            if self.save_results:
                self.save_model(
                    output_model, suffix='wfscmb', output_file=outfile, format=False
                )

        # Short-circuit auto-save of returned model if run from strun, as it is
        # already done above.  Ideally we would use self.output_use_model,
        # self.suffix and output_model.meta.filename, but self.save_model(format=False)
        # also needs to be there, and that is not stored as a class/instance attr.
        self.save_results = False
        # Return the output so it can be tested.  Assumes there is only one product.
        return output_model
