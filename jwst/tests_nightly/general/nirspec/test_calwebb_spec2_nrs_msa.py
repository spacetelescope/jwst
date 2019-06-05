"""Test calwebb_spec2 for NIRSpec MSA"""
import os.path as op

import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations import load_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step


@pytest.mark.bigdata
class TestSpec2NRSMSA(BaseJWSTTest):
    """Test various aspects of calibrating NIRSpec MSA mode"""
    input_loc = 'nirspec'
    ref_loc = ['test_datasets', 'msa', 'simulated-3nod', 'truth']
    test_dir = ['test_datasets', 'msa', 'simulated-3nod']

    def test_msa_missing(self, caplog):
        """Test MSA missing failure"""
        input_file = self.get_data(
            *self.test_dir, 'level2a_twoslit', 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
        )

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            input_file
        ]

        with pytest.raises(Exception):
            Step.from_cmdline(args)

        assert 'Missing MSA meta (MSAMETFL) file' in caplog.text

    def test_msa_missing_nofail(self, caplog):
        """Test MSA missing failure"""
        input_file = self.get_data(
            *self.test_dir, 'level2a_twoslit', 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
        )

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            input_file,
            '--fail_on_exception=false'
        ]

        Step.from_cmdline(args)

        assert 'Missing MSA meta (MSAMETFL) file' in caplog.text

    def test_msa_missing_skip(self, caplog):
        """Test MSA missing failure"""
        input_file = self.get_data(
            *self.test_dir, 'level2a_twoslit', 'F170LP-G235M_MOS_observation-6-c0e0_001_DN_NRS1_mod.fits'
        )

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            input_file,
            '--steps.assign_wcs.skip=true'
        ]

        Step.from_cmdline(args)

        assert 'Aborting remaining processing for this exposure.' in caplog.text

    def test_run_msaflagging(self, caplog):
        """Test msa flagging operation"""

        # Retrieve the data.
        collect_pipeline_cfgs('cfgs')
        self.get_data(
            *self.test_dir, 'jw95065006001_0_msa_twoslit.fits'
        )
        asn_path = self.get_data(
            *self.test_dir, 'mos_udf_g235m_twoslit_spec2_asn.json'
        )
        with open(asn_path) as fp:
            asn = load_asn(fp)
        for product in asn['products']:
            for member in product['members']:
                self.get_data(
                    *self.test_dir, 'level2a_twoslit', member['expname']
                )

        # Run step.
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            asn_path,
            '--steps.msa_flagging.skip=false'
        ]
        Step.from_cmdline(args)

        # Test.
        assert 'Step msa_flagging running with args' in caplog.text
        assert 'Step msa_flagging done' in caplog.text

        for product in asn['products']:
            prod_name = product['name'] + '_cal.fits'
            assert op.isfile(prod_name)
