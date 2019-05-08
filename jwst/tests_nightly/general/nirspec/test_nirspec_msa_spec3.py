"""Test calwebb_spec3 against NIRSpec MOS science (MSA)"""
from pathlib import Path
import pytest

from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step

from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn


@pytest.mark.bigdata
class TestSpec3Pipeline(BaseJWSTTest):
    """Tests for Spec3Pipeline"""

    input_loc = 'nirspec'
    ref_loc = ['test_datasets', 'msa', 'sdp_jw95175', 'truth']
    test_dir = ['test_datasets', 'msa', 'sdp_jw95175']

    def test_nrs_msa_spec3(self):
        """
        Regression test of calwebb_spec3 pipeline performed on
        NIRSpec MSA data
        """
        cfg_dir = './cfgs'
        collect_pipeline_cfgs(cfg_dir)
        asn_file = self.get_data(*self.test_dir,
                                 'single_asn.json')

        for file in raw_from_asn(asn_file):
            self.get_data(*self.test_dir, file)

        args = [
            str(Path(cfg_dir) / 'calwebb_spec3.cfg'),
            asn_file
        ]

        Step.from_cmdline(args)

        # Compare results
        truths = self.data_glob(*self.ref_loc, glob='*.fits')
        outputs = [
            (Path(output_file).name, ) * 2
            for output_file in truths
        ]
        self.compare_outputs(outputs)
