import pytest
from jwst.tests.base_test import BaseJWSTTest

from jwst.pipeline.calwebb_image3 import Image3Pipeline


@pytest.mark.bigdata
class TestImage3Pipeline1(BaseJWSTTest):
    """Regression test definitions for CALIMAGE3 pipeline.

    Regression test of calwebb_image3 pipeline on NIRCam
    simulated long-wave data.
    """
    input_loc = 'nircam'
    ref_loc = ['test_calimage3', 'truth']

    def test_image3_pipeline1(self, _jail):

        asn_name = "mosaic_long_asn.json"
        asn_file = self.get_data('test_calimage3', asn_name)
        for file in self.raw_from_asn(asn_file):
            self.get_data('test_calimage3', file)

        pipe = Image3Pipeline()
        pipe.tweakreg.skip = True
        pipe.skymethod = 'global+match'
        pipe.skymatch.match_down = True
        pipe.skymatch.subtract = False
        pipe.skymatch.skystat = 'mode'
        pipe.skymatch.nclip = 5
        pipe.skymatch.lsigma = 4.0
        pipe.skymatch.usigma = 4.0
        pipe.skymatch.binwidth = 0.1
        pipe.outlier_detection.weight_type = 'exptime'
        pipe.outlier_detection.pixfrac = 1.0
        pipe.outlier_detection.kernel = 'square'
        pipe.outlier_detection.fillval = 'INDEF'
        pipe.outlier_detection.nlow = 0
        pipe.outlier_detection.nhigh = 0
        pipe.outlier_detection.maskpt = 0.7
        pipe.outlier_detection.grow = 1
        pipe.outlier_detection.snr = '4.0 3.0'
        pipe.outlier_detection.scale = '0.5 0.4'
        pipe.outlier_detection.backg = 0.0
        pipe.outlier_detection.save_intermediate_results = False
        pipe.outlier_detection.resample_data = True
        pipe.outlier_detection.good_bits = 4
        pipe.resample.single = False
        pipe.resample.weight_type = 'exptime'
        pipe.resample.pixfrac = 1.0
        pipe.resample.kernel = 'square'
        pipe.resample.fillval = 'INDEF'
        pipe.resample.good_bits = 4
        pipe.resample.blendheaders = True
        pipe.source_catalog.kernel_fwhm = 3.
        pipe.source_catalog.kernel_xsize = 5.
        pipe.source_catalog.kernel_ysize = 5.
        pipe.source_catalog.snr_threshold = 3.
        pipe.source_catalog.npixels = 50
        pipe.source_catalog.deblend = False

        pipe.run(asn_file)

        outputs = [(# Compare level-2c crf product
                             'nrca5_47Tuc_subpix_dither1_newpos_a3001_crf.fits',
                             'nrca5_47Tuc_subpix_dither1_newpos_cal-a3001_ref.fits'),
                   {'files':(# Compare i2d product
                             'mosaic_long_i2d.fits',
                             'mosaic_long_i2d_ref.fits'),
                    'pars': {'ignore_hdus':self.ignore_hdus+['HDRTAB'],
                              'rtol': 0.0001}
                   },
                   {'files':(# Compare the HDRTAB in the i2d product
                             'mosaic_long_i2d.fits[HDRTAB]',
                             'mosaic_long_i2d_ref.fits[HDRTAB]'),
                    'pars': {'ignore_keywords':
                             self.ignore_keywords+['NAXIS1', 'TFORM*'],
                             'rtol': 0.0001}
                   }
                  ]
        self.compare_outputs(outputs)

    def test_image3_pipeline2(self, _jail):
        """Regression test definitions for CALIMAGE3 pipeline.

        Regression test of calwebb_image3 pipeline on NIRCam
        simulated long-wave data with a 6-point dither.
        """
        asn_file = self.get_data(self.test_dir,
                                 "jw10002-o001_20171116t191235_image3_002_asn.json")
        for file in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        pipe = Image3Pipeline()
        pipe.tweakreg.save_catalogs = False
        pipe.tweakreg.catalog_format = 'ecsv'
        pipe.tweakreg.kernel_fwhm = 2.
        pipe.tweakreg.snr_threshold = 5.
        pipe.tweakreg.enforce_user_order = True
        pipe.tweakreg.expand_refcat = False
        pipe.tweakreg.minobj = 15
        pipe.tweakreg.searchrad = 10.0
        pipe.tweakreg.use2dhist = True
        pipe.tweakreg.separation = 0.5
        pipe.tweakreg.tolerance = 1.0
        pipe.tweakreg.xoffset = 0.0
        pipe.tweakreg.yoffset = 0.0
        pipe.tweakreg.fitgeometry = 'rscale'
        pipe.tweakreg.nclip = 3
        pipe.tweakreg.sigma = 3.0
        pipe.skymethod = 'global+match'
        pipe.skymatch.match_down = True
        pipe.skymatch.subtract = False
        pipe.skymatch.skystat = 'mode'
        pipe.skymatch.nclip = 5
        pipe.skymatch.lsigma = 4.0
        pipe.skymatch.usigma = 4.0
        pipe.skymatch.binwidth = 0.1
        pipe.outlier_detection.weight_type = 'exptime'
        pipe.outlier_detection.pixfrac = 1.0
        pipe.outlier_detection.kernel = 'square'
        pipe.outlier_detection.fillval = 'INDEF'
        pipe.outlier_detection.nlow = 0
        pipe.outlier_detection.nhigh = 0
        pipe.outlier_detection.maskpt = 0.7
        pipe.outlier_detection.grow = 1
        pipe.outlier_detection.snr = '4.0 3.0'
        pipe.outlier_detection.scale = '0.5 0.4'
        pipe.outlier_detection.backg = 0.0
        pipe.outlier_detection.save_intermediate_results = False
        pipe.outlier_detection.resample_data = True
        pipe.outlier_detection.good_bits = 4
        pipe.resample.single = False
        pipe.resample.weight_type = 'exptime'
        pipe.resample.pixfrac = 1.0
        pipe.resample.kernel = 'square'
        pipe.resample.fillval = 'INDEF'
        pipe.resample.good_bits = 4
        pipe.resample.blendheaders = True
        pipe.resample.suffix = 'i2d'
        pipe.source_catalog.kernel_fwhm = 3.
        pipe.source_catalog.kernel_xsize = 5.
        pipe.source_catalog.kernel_ysize = 5.
        pipe.source_catalog.snr_threshold = 3.
        pipe.source_catalog.npixels = 50
        pipe.source_catalog.deblend = False

        pipe.run(asn_file)

        outputs = [('jw10002001001_01101_00004_nrcblong_o001_crf.fits',
                    'jw10002001001_01101_00004_nrcblong_o001_crf_ref.fits'),
                   {'files':('jw10002-o001_t002_nircam_f444w_i2d.fits',
                    'jw10002-o001_t002_nircam_f444w_i2d_ref.fits'),
                    'pars':{'ignore_hdus':self.ignore_hdus+['HDRTAB'],
                            'rtol':0.0001},
                   {'files':('jw10002-o001_t002_nircam_f444w_i2d.fits[hdrtab]',
                    'jw10002-o001_t002_nircam_f444w_i2d_ref.fits[hdrtab]'),
                    'pars':{'ignore_keywords':self.ignore_keywords+['NAXIS1','TFORM*'],
                            'rtol':0.0001}
                  ]
        self.compare_outputs(outputs)
