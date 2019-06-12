"""Test calwebb_spec3 against NIRSpec Fixed-slit science (FSS)"""
from glob import glob
from os import path
import pytest

from jwst.associations import load_asn
from jwst.pipeline import Spec3Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn


@pytest.mark.bigdata
class TestSpec3Pipeline(BaseJWSTTest):
    input_loc = 'nirspec'

    def test_save_source_only(self):
        """Test saving the source-based files only"""
        datapath = ['test_datasets', 'fss', '93045', 'level2b']

        asn_file = self.get_data(*datapath,
                                    'jw93045-o010_20180725t035735_spec3_001_asn.json')
        for file in raw_from_asn(asn_file):
            self.get_data(*datapath, file)

        pipe = Spec3Pipeline()
        pipe.mrs_imatch.skip = True
        pipe.outlier_detection.skip = True
        pipe.resample_spec.skip = True
        pipe.cube_build.skip = True
        pipe.extract_1d.skip = True

        pipe.run(asn_file)

        # Check resulting product
        with open(asn_file) as fh:
            asn = load_asn(fh)
        base_name = asn['products'][0]['name']
        product_name = base_name.format(source_id='s00000') + '_cal.fits'
        output_files = glob('*')

        if product_name in output_files:
            output_files.remove(product_name)
        else:
            assert False


    @pytest.mark.xfail(
        reason='Dataset fails at outlier_detection',
        run=False
    )
    def test_nrs_fs_spec3(self):
        """
        Regression test of calwebb_spec3 pipeline performed on
        NIRSpec fixed-slit data.
        """
        cfg_dir = './cfgs'
        collect_pipeline_cfgs(cfg_dir)
        datapath = ['test_datasets', 'fss', '93045', 'level2b']
        asn_file = self.get_data(*datapath,
                                 'jw93045-o010_20180725t035735_spec3_001_asn.json')

        for file in raw_from_asn(asn_file):
            self.get_data(*datapath, file)

        args = [
            path.join(cfg_dir, 'calwebb_spec3.cfg'),
            asn_file
        ]

        Step.from_cmdline(args)

        # Compare results
        outputs = [('jw00023001001_01101_00001_NRS1_cal.fits',
                    'jw00023001001_01101_00001_NRS1_cal_ref.fits'),
                   ('jw00023001001_01101_00001_NRS1_s2d.fits',
                    'jw00023001001_01101_00001_NRS1_s2d_ref.fits'),
                   ('jw00023001001_01101_00001_NRS1_x1d.fits',
                    'jw00023001001_01101_00001_NRS1_x1d_ref.fits')
                  ]
        self.compare_outputs(outputs)
