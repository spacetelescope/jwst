"""Test calwebb_tso3 with NIRISS SOSS"""

from collections import defaultdict
from glob import glob
import os.path as op

import pytest

from jwst.pipeline.tests.helpers import (
    SCRIPT_DATA_PATH,
)
from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import (Step, remove_suffix)
from jwst.white_light import WhiteLightStep


@pytest.mark.bigdata
class TestTSO3NIRSOSS(BaseJWSTTest):
    """Test various aspects of calibrating NIRISS SOSS data"""
    input_loc = 'niriss'
    ref_loc = ['test_caltso3', 'truth']
    test_dir = ['test_caltso3']

    def test_run_full_noscale(self):
        """Test a full run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw87600-a3001_20170527t111213_tso3_001_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, member['expname']
                )
        input_files = glob('*')

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_tso3.cfg'),
            asn_path,
            '--scale_detection=False',
            '--steps.outlier_detection.save_intermediate_results=True',
        ]
        Step.from_cmdline(args)

        # Now test for file existence. Get the association
        with open(asn_path) as fh:
            asn = load_asn(fh)
        acid = asn['asn_id']
        product = asn['products'][0]
        product_name = product['name']
        members_by_type = defaultdict(list)
        for member in product['members']:
            expname = op.split(member['expname'])[1]
            members_by_type[member['exptype'].lower()].append(expname)

        output_files = glob('*')
        for input_file in input_files:
            output_files.remove(input_file)
        print('Created files ares: {}'.format(output_files))

        # Check Level3 products
        product_name_file = product_name + '_whtlt.ecsv'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        product_name_file = product_name + '_x1dints.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        # Check Level2 products
        for member in members_by_type['science']:
            basename, ext = op.splitext(op.split(member)[1])
            basename, separator = remove_suffix(basename)

            name = basename + separator + acid + separator + 'crfints' + ext
            assert name in output_files
            output_files.remove(name)

            name = basename + separator + acid + separator + 'median' + ext
            assert name in output_files
            output_files.remove(name)

        # If there are files left, this is an error
        assert len(output_files) == 0

    def test_run_whitelight(self):
        """Test the whiteline step"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw87600-a3001_20170527t111213_tso3_001_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, member['expname']
                )
        input_files = glob('*')

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_tso3.cfg'),
            asn_path,
            '--steps.outlier_detection.skip=true',
            '--steps.tso_photometry.skip=true',
        ]
        Step.from_cmdline(args)

        # Now test for file existence. Get the association
        with open(asn_path) as fh:
            asn = load_asn(fh)
        product = asn['products'][0]
        product_name = product['name']
        members_by_type = defaultdict(list)
        for member in product['members']:
            expname = op.split(member['expname'])[1]
            members_by_type[member['exptype'].lower()].append(expname)

        output_files = glob('*')
        for input_file in input_files:
            output_files.remove(input_file)
        print('Created files ares: {}'.format(output_files))

        # Check Level3 products
        product_name_file = product_name + '_whtlt.ecsv'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        product_name_file = product_name + '_x1dints.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        # If there are files left, this is an error
        assert len(output_files) == 0

    def test_whitelight_output(self):
        """Test default file output from white_light_step"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        input_file = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits'
        input_path = self.get_data(
            *self.test_dir, input_file)
        input_files = glob('*')

        # Run the step.
        args = [
            op.join('cfgs', 'white_light.cfg'),
            input_path
        ]
        Step.from_cmdline(args)

        # Test.
        output_files = glob('*')
        for input_file in input_files:
            output_files.remove(input_file)
        print('Created files ares: {}'.format(output_files))

        basename, ext = op.splitext(input_file)
        basename, separator = remove_suffix(basename)
        output_name = basename + '_whtlt.ecsv'
        assert output_name in output_files
        output_files.remove(output_name)

        assert len(output_files) == 0

    def test_whitelight_nooutput(self):
        """Test for no output from white_light_step"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        input_file = 'jw87600-a3001_t1_niriss_clear-gr700xd_x1dints.fits'
        input_path = self.get_data(
            *self.test_dir, input_file)
        input_files = glob('*')

        # Run the step.
        WhiteLightStep.call(
            input_path,
            config_file=op.join('cfgs', 'white_light.cfg')
        )

        # Test.
        output_files = glob('*')
        for input_file in input_files:
            output_files.remove(input_file)
        print('Created files ares: {}'.format(output_files))

        assert len(output_files) == 0
