import pytest

from jwst.tests.base_classes import BaseJWSTTest, raw_from_asn
from jwst.pipeline.collect_pipeline_cfgs import collect_pipeline_cfgs
from jwst.stpipe import Step


@pytest.mark.bigdata
class TestImage3Pipeline1(BaseJWSTTest):
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data.
    """
    input_loc = 'nircam'
    ref_loc = ['test_calimage3', 'truth']
    test_dir = 'test_calimage3'

    def test_image3_pipeline1(self):

        asn_name = "mosaic_long_asn.json"
        asn_file = self.get_data('test_calimage3', asn_name)
        for file in raw_from_asn(asn_file):
            self.get_data('test_calimage3', file)

        collect_pipeline_cfgs('config')

        args = [
            'config/calwebb_image3.cfg',
            asn_file,
            '--steps.tweakreg.skip=True',
        ]

        Step.from_cmdline(args)

        outputs = [(# Compare level-2c crf product
                             'nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits',
                             'nrca5_47Tuc_subpix_dither1_newpos_cal-a3001_ref.fits'),
                   {'files':(# Compare i2d product
                             'mosaic_long_i2d.fits',
                             'mosaic_long_i2d_ref.fits'),
                    'pars': {'ignore_keywords':
                             self.ignore_keywords+['NAXIS1', 'TFORM*'],
                             'ignore_fields':self.ignore_keywords,
                             'rtol': 0.0001}
                   }
                  ]
        self.compare_outputs(outputs)

    def test_image3_pipeline2(self):
        """Regression test definitions for CALIMAGE3 pipeline.

        Regression test of calwebb_image3 pipeline on NIRCam
        simulated long-wave data with a 6-point dither.
        """
        asn_file = self.get_data(self.test_dir,
                                 "jw10002-o001_20171116t191235_image3_002_asn.json")
        for file in raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        collect_pipeline_cfgs('config')

        args = [
            'config/calwebb_image3.cfg',
            asn_file,
            '--steps.tweakreg.kernel_fwhm=2',
            '--steps.tweakreg.snr_threshold=5',
            '--steps.tweakreg.enforce_user_order=True',
            '--steps.tweakreg.searchrad=10',
            '--steps.tweakreg.fitgeometry=rscale'
        ]

        Step.from_cmdline(args)

        outputs = [('jw10002001001_01101_00004_nrcblong_o001_crf.fits',
                    'jw10002001001_01101_00004_nrcblong_o001_crf_ref.fits'),
                   {'files':('jw10002-o001_t002_nircam_f444w_i2d.fits',
                    'jw10002-o001_t002_nircam_f444w_i2d_ref.fits'),
                    'pars':{'ignore_keywords':self.ignore_keywords+['NAXIS1','TFORM*'],
                            'ignore_fields':self.ignore_keywords,
                            'rtol':0.0001}
                   }
                  ]
        self.compare_outputs(outputs)
