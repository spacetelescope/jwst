"""Test calwebb_spec3 with MIRI MRS with 4-dither pattern"""

import os.path as op

import pytest

from jwst.tests.base_classes import BaseJWSTTest
from jwst.pipeline.tests.helpers import SCRIPT_DATA_PATH

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step


@pytest.mark.bigdata
class TestSpec3MIRMRSDither(BaseJWSTTest):
    """Test various aspects of calibrating MIRI MRS dither data"""
    input_loc = 'miri'
    ref_loc = ['test_datasets', 'mrs', 'cv3_ifu_dither', 'truth']
    test_dir = ['test_datasets', 'mrs', 'cv3_ifu_dither']

    @pytest.mark.xfail(
        reason='Data needs reformatting: Input data is not a IFUImageModel',
        run=False
    )
    def test_run_cube_build_only(self):
        """Test only the extraction step.
        """

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'level2b', 'cube_build_4dither_495_asn.json'
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
            '--steps.extract_1d.skip=true',
        ]

        Step.from_cmdline(args)

        # Test
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_ch1-short_s3d.fits'
        assert op.isfile(product_name)
        product_name = product_name_base + '_ch2-short_s3d.fits'
        assert op.isfile(product_name)

    @pytest.mark.xfail(
        reason='Data needs reformatting: Input data is not a IFUImageModel',
        run=False
    )
    def test_run_extract_1d_only(self):
        """Test only the extraction step.
        """

        # Get the data.
        collect_pipeline_cfgs('cfgs')
        asn_path = self.get_data(
            *self.test_dir, 'level2b', 'cube_build_4dither_495_asn.json'
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
            op.join(SCRIPT_DATA_PATH, 'calwebb_spec3_default.cfg'),
            asn_path,
            '--steps.mrs_imatch.skip=true',
            '--steps.outlier_detection.skip=true',
            '--steps.resample_spec.skip=true',
        ]

        Step.from_cmdline(args)

        # Test.
        with open(asn_path) as fd:
            asn = load_asn(fd)
        product_name_base = asn['products'][0]['name']
        product_name = product_name_base + '_ch1-short_s3d.fits'
        assert op.isfile(product_name)
        product_name = product_name_base + '_ch2-short_s3d.fits'
        assert op.isfile(product_name)
        product_name = product_name_base + '_ch1-short_x1d.fits'
        assert op.isfile(product_name)
        product_name = product_name_base + '_ch2-short_x1d.fits'
        assert op.isfile(product_name)
