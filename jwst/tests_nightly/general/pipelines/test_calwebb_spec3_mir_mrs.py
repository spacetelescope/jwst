"""Test calwebb_spec3 with MIRI MRS"""

from collections import defaultdict
from glob import glob
import os.path as op

import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.pipeline.tests.helpers import (
    SCRIPT_DATA_PATH,
)

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import (Step, remove_suffix)


@pytest.mark.bigdata
class TestSpec3MIRMRS(BaseJWSTTest):
    """Test various aspects of calibrating MIRI MRS data"""
    input_loc = 'miri'
    ref_loc = ['test_datasets', 'mrs', 'simulated', 'truth']
    test_dir = ['test_datasets', 'mrs', 'simulated']

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_nothing(self):
        """Run no steps. There should be no output."""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true'
        ]
        Step.from_cmdline(args)

        # Test.
        assert len(glob('*')) == 0

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_extract_1d_only(self):
        """Test only the extraction step.
        """

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name_glob = product_name_base + '_ch[34]-long_s3d.fits'
        assert len(glob(product_name_glob)) == 2
        product_name_glob = product_name_base + '_ch[34]-long_x1d.fits'
        assert len(glob(product_name_glob)) == 2

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_resample_only(self):
        """Test resample step only."""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name_glob = product_name_base + '_ch[34]-long_s3d.fits'
        assert len(glob(product_name_glob)) == 2

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_mrs_imatch_only(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
            '--steps.mrs_imatch.save_results=true',
        ]
        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_mrs_imatch.fits'
        assert op.isfile(product_name)

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_full(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_channel_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
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

        output_files = [
            op.split(result_path)[1]
            for result_path in
            glob('*')
        ]
        print('Created files ares: {}'.format(output_files))

        # Check Level3 products
        level3_suffixes = [
            '_ch1-long_s3d.fits',
            '_ch2-long_s3d.fits',
            '_ch3-long_s3d.fits',
            '_ch4-long_s3d.fits',
            '_ch1-long_x1d.fits',
            '_ch2-long_x1d.fits',
            '_ch3-long_x1d.fits',
            '_ch4-long_x1d.fits',
        ]
        for level3_suffix in level3_suffixes:
            product_name_file = product_name + level3_suffix
            assert product_name_file in output_files
            output_files.remove(product_name_file)

        # Check Level2 products
        for member in members_by_type['science']:
            basename, ext = op.splitext(op.split(member)[1])
            basename, separator = remove_suffix(basename)

            name = basename + separator + acid + separator + 'crf' + ext
            assert name in output_files
            output_files.remove(name)

        # If there are files left, this is an error
        assert len(output_files) == 0

    @pytest.mark.xfail(
        reason="Data needs reformatting: TaggedDict 'TaggedDict' object has no attribute '_naxes'",
        run=False,
    )
    def test_run_outlier_only(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'single_channel_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Now test for file existence. Get the association
        with open(asn_path) as fh:
            asn = load_asn(fh)
        acid = asn['asn_id']
        product = asn['products'][0]
        members_by_type = defaultdict(list)
        for member in product['members']:
            expname = op.split(member['expname'])[1]
            members_by_type[member['exptype'].lower()].append(expname)

        output_files = [
            op.split(result_path)[1]
            for result_path in
            glob('*')
        ]
        print('Created files ares: {}'.format(output_files))

        # Check Level2 products
        for member in members_by_type['science']:
            basename, ext = op.splitext(op.split(member)[1])
            basename, separator = remove_suffix(basename)

            name = basename + separator + acid + separator + 'crf' + ext
            assert name in output_files
            output_files.remove(name)

        # If there are files left, this is an error
        assert len(output_files) == 0
