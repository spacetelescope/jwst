"""Test wfs_combine"""

from glob import glob
import os.path as op

import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.associations.lib.rules_level3_base import format_product
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step


@pytest.mark.bigdata
class TestWFSImage3Pipeline(BaseJWSTTest):
    input_loc = 'nircam'
    ref_loc = ['test_wfs_combine', 'truth']
    test_dir = 'test_wfs_combine'

    def test_asn_naming(self):
        """Test a full run"""

        # Get the data
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            self.test_dir, 'wfs_3sets_asn.json'
        )
        with open(asn_path) as fh:
            asn = load_asn(fh)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    self.test_dir, member['expname']
                )
        input_files = glob('*')

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_wfs-image3.cfg'),
            asn_path
        ]
        Step.from_cmdline(args)

        # Test.
        output_files = glob('*')
        for input_file in input_files:
            output_files.remove(input_file)
        print('output_files = {}'.format(output_files))

        for product in asn['products']:
            prod_name = product['name']
            prod_name = format_product(prod_name, suffix='wfscmb')
            prod_name += '.fits'
            assert prod_name in output_files
            output_files.remove(prod_name)

        # There should be no more files
        assert len(output_files) == 0
