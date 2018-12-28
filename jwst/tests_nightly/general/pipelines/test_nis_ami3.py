"""Test calwebb_ami3 with NIR"""

from collections import defaultdict
from glob import glob
import os.path as op

import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import (Step, remove_suffix)


@pytest.mark.bigdata
class TestAMI3NIR(BaseJWSTTest):
    """Test various aspects of calibrating NIRISS AMI data"""
    input_loc = 'niriss'
    ref_loc = ['test_ami_pipeline', 'truth']
    test_dir = ['test_ami_pipeline']

    def test_run_full(self):
        """Test a full run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'test_lg1_asn.json'
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
            op.join('cfgs', 'calwebb_ami3.cfg'),
            asn_path,
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
        product_name_file = product_name + '_amiavg.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        product_name_file = product_name + '_psf-amiavg.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        product_name_file = product_name + '_aminorm.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        # Check Level2 products
        for member in members_by_type['psf']:
            name, ext = op.splitext(op.split(member)[1])
            name, separator = remove_suffix(name)
            name = name + separator + acid + separator + 'ami' + ext
            assert name in output_files
            output_files.remove(name)

        for member in members_by_type['science']:
            name, ext = op.splitext(op.split(member)[1])
            name, separator = remove_suffix(name)
            name = name + separator + acid + separator + 'ami' + ext
            assert name in output_files
            output_files.remove(name)

        # If there are files left, this is an error
        assert len(output_files) == 0

    def test_run_single(self):
        """Run a single exposure"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'test_lg1_single_asn.json'
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
            op.join('cfgs', 'calwebb_ami3.cfg'),
            asn_path,
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
        product_name_file = product_name + '_amiavg.fits'
        assert product_name_file in output_files
        output_files.remove(product_name_file)

        # Check Level2 products
        for member in members_by_type['psf']:
            name, ext = op.splitext(op.split(member)[1])
            name, separator = remove_suffix(name)
            name = name + separator + acid + separator + 'ami' + ext
            assert name in output_files
            output_files.remove(name)

        for member in members_by_type['science']:
            name, ext = op.splitext(op.split(member)[1])
            name, separator = remove_suffix(name)
            name = name + separator + acid + separator + 'ami' + ext
            assert name in output_files
            output_files.remove(name)

        # If there are files left, this is an error
        assert len(output_files) == 0
