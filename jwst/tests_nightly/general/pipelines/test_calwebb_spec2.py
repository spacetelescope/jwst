"""Test calwebb_spec2"""

from glob import glob
import os.path as op
import pytest

from jwst.tests.base_classes import BaseJWSTTest

from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.datamodels import DataModel
from jwst.datamodels import open as dm_open
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe.step import Step

EXPFILE = 'jw00035001001_01101_00001_MIRIMAGE_rate.fits'
CALFILE = EXPFILE.replace('_rate', '_cal')
BSUBFILE = EXPFILE.replace('_rate', '_bsub')
EXTRACT1DFILE = EXPFILE.replace('_rate', '_x1d')


@pytest.mark.bigdata
class TestCalwebbSpec2(BaseJWSTTest):
    """Test various aspects of the calwebb_spec2 pipeline"""
    input_loc = 'miri'
    ref_loc = ['test_spec2pipeline', 'truth']
    test_dir = 'test_spec2pipeline'

    def test_result_return(self):
        """Ensure that a result is returned programmatically"""
        collect_pipeline_cfgs('cfgs')
        exppath = self.get_data(self.test_dir, EXPFILE)

        results = Spec2Pipeline.call(exppath, config_file=op.join('cfgs', 'calwebb_spec2.cfg'))
        assert isinstance(results[0], DataModel)


    def test_full_run(self):
        """Make a full run with the default configuration"""
        collect_pipeline_cfgs('cfgs')
        exppath = self.get_data(self.test_dir, EXPFILE)

        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            exppath
        ]
        Step.from_cmdline(args)

        assert op.isfile(CALFILE)
        assert op.isfile(EXTRACT1DFILE)


    def test_asn_with_bkg(self):
        """Test association input with background subtraction"""
        # Setup the association
        lv2_meta = {
            'program': 'test',
            'target': 'test',
            'asn_pool': 'test',
        }
        exppath = self.get_data(self.test_dir, EXPFILE)
        asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
        asn['products'][0]['members'].append({
            'expname': exppath, 'exptype': 'BACKGROUND'
        })
        asn_path, serialized = asn.dump()
        with open(asn_path, 'w') as fp:
            fp.write(serialized)

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            asn_path,
            '--steps.bkg_subtract.save_results=true'
        ]

        Step.from_cmdline(args)

        output_files = glob('*')
        print('Created files are: {}'.format(output_files))

        assert op.isfile(CALFILE)
        assert op.isfile(BSUBFILE)


    def test_asn_with_bkg_bsub(self):
        """Test background saving using the bsub flag"""
        lv2_meta = {
            'program': 'test',
            'target': 'test',
            'asn_pool': 'test',
        }
        exppath = self.get_data(self.test_dir, EXPFILE)
        asn = asn_from_list([exppath], rule=DMSLevel2bBase, meta=lv2_meta)
        asn['products'][0]['members'].append({
            'expname': exppath, 'exptype': 'BACKGROUND'
        })
        asn_file, serialized = asn.dump()
        with open(asn_file, 'w') as fp:
            fp.write(serialized)

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            asn_file,
            '--save_bsub=true'
        ]

        Step.from_cmdline(args)

        assert op.isfile(CALFILE)
        assert op.isfile(BSUBFILE)


    def test_asn(self):
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

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            asn_file,
        ]

        Step.from_cmdline(args)

        assert op.isfile(CALFILE)


    def test_asn_multiple_products(self):
        """Test an association with multiple products"""
        lv2_meta = {
            'program': 'test',
            'target': 'test',
            'asn_pool': 'test',
        }
        exppath = self.get_data(self.test_dir, EXPFILE)
        asn = asn_from_list([exppath, exppath], rule=DMSLevel2bBase, meta=lv2_meta)
        asn['products'][0]['name'] = 'product1'
        asn['products'][1]['name'] = 'product2'
        asn_file, serialized = asn.dump()
        with open(asn_file, 'w') as fp:
            fp.write(serialized)

        collect_pipeline_cfgs('cfgs')
        args = [
            op.join('cfgs', 'calwebb_spec2.cfg'),
            asn_file,
        ]

        Step.from_cmdline(args)

        assert op.isfile('product1_cal.fits')
        assert op.isfile('product2_cal.fits')


    def test_datamodel(self):
        """Test using a DataModel as input"""
        model = dm_open(self.get_data(self.test_dir, EXPFILE))
        collect_pipeline_cfgs('cfgs')
        Spec2Pipeline.call(model, config_file=op.join('cfgs', 'calwebb_spec2.cfg'))
        assert op.isfile(CALFILE)


    def test_file(self):
        """Test using a file as input"""
        exppath = self.get_data(self.test_dir, EXPFILE)
        collect_pipeline_cfgs('cfgs')
        Spec2Pipeline.call(exppath, config_file=op.join('cfgs', 'calwebb_spec2.cfg'))
        assert op.isfile(CALFILE)
