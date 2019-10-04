import pytest

from jwst.pipeline import Spec2Pipeline

from jwst.tests.base_classes import BaseJWSTTest
from jwst.tests.base_classes import pytest_generate_tests # noqa: F401


@pytest.mark.bigdata
class TestSpec2Pipeline(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_pipelines', 'truth']
    test_dir = 'test_pipelines'

    # Specification of parameters for Spec2Pipeline tests
    params = {'test_spec2':
                # test_nrs_fs_multi_spec2_1: NIRSpec fixed-slit data
                [dict(input='jw00023001001_01101_00001_NRS1_rate.fits',
                      outputs=[('jw00023001001_01101_00001_NRS1_cal.fits',
                                  'jw00023001001_01101_00001_NRS1_cal_ref.fits'),
                                 ('jw00023001001_01101_00001_NRS1_s2d.fits',
                                  'jw00023001001_01101_00001_NRS1_s2d_ref.fits'),
                                 ('jw00023001001_01101_00001_NRS1_x1d.fits',
                                  'jw00023001001_01101_00001_NRS1_x1d_ref.fits')
                                 ],
                       id="nirspec_fs_multi_1"
                      ),
                # test_nrs_fs_multi_spec2_2: NIRSpec fixed-slit data
                 dict(input= 'jwtest1013001_01101_00001_NRS1_rate.fits',
                      outputs=[('jwtest1013001_01101_00001_NRS1_cal.fits',
                                'jwtest1013001_01101_00001_NRS1_cal_ref.fits'),
                               ('jwtest1013001_01101_00001_NRS1_s2d.fits',
                                'jwtest1013001_01101_00001_NRS1_s2d_ref.fits'),
                               ('jwtest1013001_01101_00001_NRS1_x1d.fits',
                                'jwtest1013001_01101_00001_NRS1_x1d_ref.fits')
                              ],
                      id="nirspec_fs_multi_2"
                     ),
                # test_nrs_fs_multi_spec2_3:
                # NIRSpec fixed-slit data using the ALLSLITS subarray and detector NRS2
                # NIRSpec fixed-slit data that uses a single-slit subarray (S200B1).
                 dict(input= 'jw84600002001_02101_00001_nrs2_rate.fits',
                      outputs=[('jw84600002001_02101_00001_nrs2_cal.fits',
                                'jw84600002001_02101_00001_nrs2_cal_ref.fits'),
                               ('jw84600002001_02101_00001_nrs2_s2d.fits',
                                'jw84600002001_02101_00001_nrs2_s2d_ref.fits'),
                               ('jw84600002001_02101_00001_nrs2_x1d.fits',
                                'jw84600002001_02101_00001_nrs2_x1d_ref.fits')
                              ],
                      id="nirspec_fs_multi_3"
                     ),
                # test_nrs_ifu_spec2: NIRSpec IFU data
                 dict(input= 'jw95175001001_02104_00001_nrs1_rate.fits',
                      outputs=[('jw95175001001_02104_00001_nrs1_cal.fits',
                                'jw95175001001_02104_00001_nrs1_cal_ref.fits'),
                               ('jw95175001001_02104_00001_nrs1_s3d.fits',
                                'jw95175001001_02104_00001_nrs1_s3d_ref.fits'),
                               ('jw95175001001_02104_00001_nrs1_x1d.fits',
                                'jw95175001001_02104_00001_nrs1_x1d_ref.fits')
                              ],
                      id = "nirspec_ifu"
                      )
                ]
            }

    def test_spec2(self, input, outputs):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec data.
        """
        input_file = self.get_data(self.test_dir, input)

        step = Spec2Pipeline()
        step.save_bsub = True
        step.save_results = True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        self.compare_outputs(outputs)
