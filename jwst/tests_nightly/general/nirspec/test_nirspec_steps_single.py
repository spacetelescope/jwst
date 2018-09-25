import pytest
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.tests.base_test import BaseJWSTTest

from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.imprint.imprint_step import ImprintStep
from jwst.ramp_fitting import RampFitStep


@pytest.mark.bigdata
class TestDetector1Pipeline(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_pipelines', 'truth']
    test_dir = 'test_pipelines'

    def test_detector1pipeline4(self):
        """

        Regression test of calwebb_detector1 pipeline performed on NIRSpec data.

        """
        input_file = self.get_data(self.test_dir,
                                   'jw84600007001_02101_00001_nrs1_uncal.fits')
        step = Detector1Pipeline()
        step.save_calibrated_ramp = True
        step.ipc.skip = True
        step.persistence.skip = True
        step.jump.rejection_threshold = 4.0
        step.ramp_fit.save_opt = False
        step.output_file = 'jw84600007001_02101_00001_nrs1_rate.fits'
        step.run(input_file)

        outputs = [('jw84600007001_02101_00001_nrs1_ramp.fits',
                    'jw84600007001_02101_00001_nrs1_ramp_ref.fits'),
                   ('jw84600007001_02101_00001_nrs1_rate.fits',
                    'jw84600007001_02101_00001_nrs1_rate_ref.fits')
                  ]
        self.compare_outputs(outputs)

@pytest.mark.bigdata
class TestNIRSpecImprint(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_imprint', 'truth']
    test_dir = 'test_imprint'

    def test_imprint_nirspec(self):
        """

        Regression test of imprint step performed on NIRSpec MSA data.

        """
        input_file = self.get_data(self.test_dir,
                                   'jw00038001001_01101_00001_NRS1_rate.fits')
        model_file = self.get_data(self.test_dir,
                                   'NRSMOS-MODEL-21_NRS1_rate.fits')

        result = ImprintStep.call(input_file, model_file, name='imprint')

        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        outputs = [(output_file,
                    'jw00038001001_01101_00001_NRS1_imprint.fits')]
        self.compare_outputs(outputs)


class TestNIRSpecRampFit(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_ramp_fit', 'truth']
    test_dir = 'test_ramp_fit'

    def test_ramp_fit_nirspec(self):
        """

        Regression test of ramp_fit step performed on NIRSpec data. This is a single
        integration dataset.

        """
        input_file = self.get_data(self.test_dir,
                                    'jw00023001001_01101_00001_NRS1_jump.fits')

        result, result_int = RampFitStep.call(input_file,
                          save_opt=True,
                          opt_name='rampfit_opt_out.fits', name='RampFit'
                          )
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        outputs = [(output_file,
                     'jw00023001001_01101_00001_NRS1_ramp_fit.fits'),
                    ('rampfit_opt_out_fitopt.fits',
                     'jw00023001001_01101_00001_NRS1_opt.fits',
                     ['primary','slope','sigslope','yint','sigyint',
                      'pedestal','weights','crmag'])
                  ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestNIRSpecWCS(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_wcs', 'nrs1-fs', 'truth']
    test_dir = ['test_wcs', 'nrs1-fs']

    def test_nirspec_nrs1_wcs(self):
        """

        Regression test of creating a WCS object and doing pixel to sky transformation.

        """
        input_file = self.get_data(*self.test_dir,
                                  'jw00023001001_01101_00001_NRS1_ramp_fit.fits')
        ref_file = self.get_data(*self.ref_loc,
                                 'jw00023001001_01101_00001_NRS1_ramp_fit_assign_wcs.fits')

        result = AssignWcsStep.call(input_file, save_results=True, suffix='assign_wcs')
        result.close()

        im = ImageModel(result.meta.filename)
        imref = ImageModel(ref_file)

        for slit in ['S200A1', 'S200A2', 'S400A1', 'S1600A1']:
            w = nirspec.nrs_wcs_set_input(im, slit)
            grid = grid_from_bounding_box(w.bounding_box)
            ra, dec, lam = w(*grid)
            wref = nirspec.nrs_wcs_set_input(imref, slit)
            raref, decref, lamref = wref(*grid)

            assert_allclose(ra, raref, equal_nan=True)
            assert_allclose(dec, decref, equal_nan=True)
            assert_allclose(lam, lamref, equal_nan=True)
