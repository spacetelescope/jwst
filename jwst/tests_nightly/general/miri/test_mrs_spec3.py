import pytest
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline

from ..resources import MIRITest

pytestmark = [
    pytest.mark.usefixtures('_jail'),
    pytest.mark.skipif(not pytest.config.getoption('bigdata'),
                       reason='requires --bigdata')
]


class TestSpec3Pipeline(MIRITest):
    ref_loc = ['mrs_calspec3']
    test_dir = 'mrs_calspec3'
    rtol = 0.000001

    def test_spec3_pipeline1(self):
        """
        Regression test of calwebb_spec3 pipeline on simulated
        MIRI MRS dithered data.
        """

        asn_file = self.get_data(self.test_dir, 'test_asn4.json')
        for file in self.raw_from_asn(asn_file):
            self.get_data(self.test_dir, file)

        step = Spec3Pipeline()
        step.save_bsub = False
        step.mrs_imatch.suffix = 'mrs_imatch'
        step.mrs_imatch.bkg_degree = 1
        step.mrs_imatch.subtract = False
        step.outlier_detection.skip = True
        step.output_use_model = True
        step.resample_spec.save_results = True
        step.resample_spec.suffix = 's2d'
        step.cube_build.save_results = True
        step.cube_build.suffix = 's3d'
        step.extract_1d.save_results = True
        step.extract_1d.suffix = 'x1d'
        step.run(asn_file)

        outputs = [(# Compare cube product 1
                    'det_image_ch1-short_s3d.fits',
                    'det_image_ch1-short_s3d_ref.fits',
                    ['primary', 'sci', 'err', 'dq', 'wmap']),
                   (# Compare cube product 2
                    'det_image_ch2-short_s3d.fits',
                    'det_image_ch2-short_s3d_ref.fits',
                    ['primary', 'sci', 'err', 'dq', 'wmap']),
                   (# Compare x1d product 1
                    'det_image_ch1-short_x1d.fits',
                    'det_image_ch1-short_x1d_ref.fits',
                    ['primary', 'extract1d']),
                   (# Compare x1d product 2
                    'det_image_ch2-short_x1d.fits',
                    'det_image_ch2-short_x1d_ref.fits',
                    ['primary', 'extract1d'])
                   ]
        self.compare_outputs(outputs)
