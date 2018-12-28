"""Test calwebb_spec3"""

from glob import glob
import os.path as op
from shutil import copy as file_copy

import pytest

from jwst.pipeline.tests.helpers import (
    SCRIPT_DATA_PATH,
)
from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step


@pytest.mark.bigdata
class TestSpec3NRSMSA(BaseJWSTTest):
    """Test various aspects of calibrating NIRSpec MSA data"""
    input_loc = 'nirspec'
    ref_loc = ['test_datasets', 'msa', 'simulated-3nod', 'truth']
    test_dir = ['test_datasets', 'msa', 'simulated-3nod']

    @pytest.mark.xfail(
        reason='Data needs reformatting: Unsupported operands for |',
        run=False,
    )
    def test_run_outlier_only(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.save_results=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        assert False

    def test_run_outlier_only_mock(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        # Check for the Source-based cal name.
        with open(asn_path) as fp:
            asn = load_asn(fp)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Check for the outlier resutls
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_crj.fits'
        assert len(glob(product_name_glob)) == 2

    @pytest.mark.xfail(
        reason="Data needs reformatting: No attribute 'to_flat_dict'",
        run=False
    )
    def test_run_resample_only(self):
        """Test resample step only."""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Check for resample results
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_s2d.fits'
        assert len(glob(product_name_glob)) == 2

    def test_run_resample_mock_only(self):
        """Test resample step only."""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Check for resample results
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_s2d.fits'
        assert len(glob(product_name_glob)) == 2

    def test_run_cube_build(self):
        """NRS MSA data is not cube data. Nothing should happen"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Check for the Source-based cal name.
        with open(asn_path) as fp:
            asn = load_asn(fp)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Assert that no cubes were built.
        cube_files = glob('*s3d*')
        assert not cube_files

    def test_run_extract_1d_only(self):
        """Test only the extraction step. Should produce nothing
        because extraction requires resampling
        """

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
        ]
        Step.from_cmdline(args)

        # Though the calibration is not run, the conversion to
        # source base has occured. Check
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Check that no other products built
        files = glob('*s3d*')
        files.extend(glob('*s2d*'))
        files.extend(glob('*x1d*'))
        assert not files

    def test_run_extract_1d_resample_mock(self):
        """Test only the extraction step. Should produce nothing
        because extraction requires resampling
        """

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.cube_build.skip=true',
        ]
        Step.from_cmdline(args)

        # Though the calibration is not run, the conversion to
        # source base has occured. Check
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_s2d.fits'
        assert len(glob(product_name_glob)) == 2

        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_x1d.fits'
        assert len(glob(product_name_glob)) == 2

    def test_run_nosteps(self):
        """Test where no steps execute"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]
        Step.from_cmdline(args)

        # Check for the Source-based cal name.
        with open(asn_path) as fp:
            asn = load_asn(fp)
        product_name_template = asn['products'][0]['name']
        product_name_glob = product_name_template.format(
            source_id='s0000[14]',
        ) + '_cal.fits'
        assert len(glob(product_name_glob)) == 2

        # Check that no other products built
        files = glob('*s3d*')
        files.extend(glob('*s2d*'))
        files.extend(glob('*x1d*'))
        assert not files

    @pytest.mark.xfail(
        reason='Data needs reformatting: Unsupported operands for |',
        run=False,
    )
    def test_run_full(self):
        """Test a full run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'two_member_spec3_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2b_twoslit', member['expname']
                )

        # Run the step.
        args = [
            op.join('cfgs', 'calwebb_spec3.cfg'),
            asn_path,
        ]

        Step.from_cmdline(args)
        assert False
