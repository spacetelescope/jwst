"""Test calwebb_image2"""
import glob
import os.path as op
from pathlib import Path

import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.datamodels import DataModel
from jwst.datamodels import open as dm_open
from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step

EXPFILE = 'jw00001001001_01101_00001_mirimage_rate.fits'
CALFILE = 'jw00001001001_01101_00001_mirimage_cal.fits'


@pytest.mark.bigdata
class TestCalwebbImage2(BaseJWSTTest):
    """Test various aspects of the calwebb_image2 pipeline"""
    input_loc = 'miri'
    ref_loc = ['test_image2pipeline', 'truth']
    test_dir = 'test_image2pipeline'

    def test_result_return(self):
        """Ensure that a result is returned programmatically"""
        collect_pipeline_cfgs('cfgs')
        exppath = self.get_data(self.test_dir, EXPFILE)

        results = Image2Pipeline.call(exppath, config_file=op.join('cfgs', 'calwebb_image2.cfg'))
        assert isinstance(results[0], DataModel)


    def test_no_cfg(self):
        """What happens when the pipeline is run without a config"""
        args = [
            'jwst.pipeline.calwebb_image2.Image2Pipeline',
            self.get_data(self.test_dir, EXPFILE)
        ]
        Step.from_cmdline(args)

        output_files = glob.glob('*')

        assert CALFILE in output_files


    def test_asn(self):
        """Test running from an association"""
        collect_pipeline_cfgs('cfgs')

        lv2_meta = {
            'program': 'test',
            'target': 'test',
            'asn_pool': 'test',
        }
        exppath = self.get_data(self.test_dir, EXPFILE)
        asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
        asn_file, serialized = asn.dump()
        with open(asn_file, 'w') as fp:
            fp.write(serialized)

        args = [
            op.join('cfgs', 'calwebb_image2.cfg'),
            asn_file,
        ]
        Step.from_cmdline(args)

        output_files = glob.glob('*')

        assert CALFILE in output_files


    def test_datamodel(self):
        """Run with a datamodel as input"""
        collect_pipeline_cfgs('cfgs')

        model = dm_open(self.get_data(self.test_dir, EXPFILE))
        Image2Pipeline.call(
            model,
            config_file=op.join('cfgs', 'calwebb_image2.cfg'))

        output_files = glob.glob('*')

        assert CALFILE in output_files


    def test_file(self):
        """Test running from an input file"""
        collect_pipeline_cfgs('cfgs')
        Image2Pipeline.call(
            self.get_data(self.test_dir, EXPFILE),
            config_file=op.join('cfgs', 'calwebb_image2.cfg'))

        output_files = glob.glob('*')

        assert CALFILE in output_files


    def test_file_outputdir(self):
        """Test putting results in another folder"""
        outfile = Path('junk.fits')
        outdir = Path('junk')
        outdir.mkdir()

        cfgs = Path('cfgs')
        collect_pipeline_cfgs(cfgs)

        args = [
            str(cfgs / 'calwebb_image2.cfg'),
            self.get_data(self.test_dir, EXPFILE),
            '--output_file=' + str(outfile),
            '--output_dir=' + str(outdir),
        ]

        Step.from_cmdline(args)

        result_path = outdir / (outfile.stem + '_cal' + outfile.suffix)
        assert result_path.exists()
