"""Test calwebb_spec3 with MIRI LRS-Fixedslit"""

from glob import glob
import os.path as op
import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step

from jwst.pipeline.tests.helpers import SCRIPT_DATA_PATH


@pytest.mark.bigdata
class TestSpec3Pipeline(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_datasets', 'lrs_fixedslit', 'truth']
    test_dir = ['test_datasets', 'lrs_fixedslit']

    @pytest.mark.xfail(
        reason='Data needs reformatting: no bounding box',
        run=False,
    )
    def test_run_outlier_only(self):
        """Test a basic run"""

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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

        # Test. 
        assert False

    def test_run_resample_mock_only(self):
        """Test resample step only."""

        # Get the data.
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]

        Step.from_cmdline(args)

        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_s2d.fits'
        assert op.isfile(product_name)

    @pytest.mark.xfail(
        reason='Data needs reformatting: bad forward transform',
        run=False
    )
    def test_run_extract_1d_resample_mock_only(self):
        """Test resample step only."""

        # Get the data.
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_mock.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.cube_build.skip=true',
        ]

        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_s2d.fits'
        assert op.isfile(product_name)
        product_name = product_name_base + '_x1d.fits'
        assert op.isfile(product_name)

    @pytest.mark.xfail(
        reason='Data needs reformatting: no bounding box',
        run=False,
    )
    def test_run_resample_only(self):
        """Test resample step only."""
        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
            '--steps.cube_build.skip=true',
            '--steps.extract_1d.skip=true',
        ]

        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_s2d.fits'
        assert op.isfile(product_name)

    def test_run_extract_1d_only(self):
        """Test only the extraction step. Should produce nothing
        because extraction requires resampling
        """
        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
        ]

        Step.from_cmdline(args)

        # Test.
        # Check that no other products built
        files = glob('*s3d*')
        files.extend(glob('*s2d*'))
        files.extend(glob('*x1d*'))
        assert not files

    def test_run_nosteps(self):
        """Test where no steps execute"""
        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
            '--steps.extract_1d.skip=true',
        ]

        Step.from_cmdline(args)

        # Check that no other products built
        files = glob('*s3d*')
        files.extend(glob('*s2d*'))
        files.extend(glob('*x1d*'))
        assert not files

    @pytest.mark.xfail(
        reason='None of the steps operate',
        run=False,
    )
    def test_run_mir_lrs_fixedslit(self):
        """Test a full run"""
        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'jw80600-a3002_20170227t160430_spec3_001_asn.json'
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
        assert False
