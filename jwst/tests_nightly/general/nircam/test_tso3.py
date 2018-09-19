from jwst.pipeline.calwebb_tso3 import Tso3Pipeline

from ..resources import NIRCamTest

class TestTso3Pipeline(NIRCamTest):
    ref_loc = ['test_caltso3']
    test_dir = 'test_caltso3'

    def test_tso3_pipeline_nrc1(self):
        """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

        Default imaging mode outlier_detection will be tested here.
        """
        asn_file = self.get_data(self.test_dir,
                                "jw93065-a3001_20170511t111213_tso3_001_asn.json")
        for file in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        step = Tso3Pipeline()
        step.scale_detection = False
        step.outlier_detection.weight_type = 'exptime'
        step.outlier_detection.pixfrac = 1.0
        step.outlier_detection.kernel = 'square'
        step.outlier_detection.fillval = 'INDEF'
        step.outlier_detection.nlow = 0
        step.outlier_detection.nhigh = 0
        step.outlier_detection.maskpt = 0.7
        step.outlier_detection.grow = 1
        step.outlier_detection.snr = '4.0 3.0'
        step.outlier_detection.scale = '0.5 0.4'
        step.outlier_detection.backg = 0.0
        step.outlier_detection.save_intermediate_results = False
        step.outlier_detection.resample_data = False
        step.outlier_detection.good_bits = 4
        step.extract_1d.smoothing_length = 0
        step.extract_1d.bkg_order = 0

        step.run(asn_file)

        outputs = [('jw93065002001_02101_00001_nrca1_a3001_crfints.fits',
                    'jw93065002001_02101_00001_nrca1_a3001_crfints_ref.fits',
                    ['primary', 'sci', 'dq', 'err'])
                   ]
        self.compare_outputs(outputs)


    def test_tso3_pipeline_nrc2(self):
        """Regression test of calwebb_tso3 pipeline on NIRCam simulated data.

        Scaled imaging mode outlier_detection will be tested here.
        """
        asn_file = self.get_data(self.test_dir,
                                 "jw93065-a3002_20170511t111213_tso3_001_asn.json")
        for file in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        step = Tso3Pipeline()
        step.scale_detection = True
        step.outlier_detection.weight_type = 'exptime'
        step.outlier_detection.pixfrac = 1.0
        step.outlier_detection.kernel = 'square'
        step.outlier_detection.fillval = 'INDEF'
        step.outlier_detection.nlow = 0
        step.outlier_detection.nhigh = 0
        step.outlier_detection.maskpt = 0.7
        step.outlier_detection.grow = 1
        step.outlier_detection.snr = '4.0 3.0'
        step.outlier_detection.scale = '0.5 0.4'
        step.outlier_detection.backg = 0.0
        step.outlier_detection.save_intermediate_results = False
        step.outlier_detection.resample_data = False
        step.outlier_detection.good_bits = 4
        step.extract_1d.smoothing_length = 0
        step.extract_1d.bkg_order = 0

        step.run(asn_file)
        outputs = [('jw93065002002_02101_00001_nrca1_a3002_crfints.fits',
                    'jw93065002002_02101_00001_nrca1_a3002_crfints_ref.fits',
                    ['primary', 'sci', 'dq', 'err'])
                  ]
        self.compare_outputs(outputs)
