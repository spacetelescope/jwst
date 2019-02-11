import pytest
from astropy.convolution import convolve, Box2DKernel
from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box

from jwst.tests.base_classes import BaseJWSTTest

from jwst.assign_wcs import AssignWcsStep, nirspec

from jwst.datamodels import ImageModel
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline
from jwst.imprint import ImprintStep
from jwst.ramp_fitting import RampFitStep
from jwst.extract_1d import Extract1dStep
#from jwst.master_background_step import MasterBackgroundStep
from jwst import datamodels


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


@pytest.mark.bigdata
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


@pytest.mark.bigdata
class TestNIRISSSpec2(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_pipelines', 'truth']
    test_dir = 'test_pipelines'

    def test_nrs_fs_single_spec2(self):
        """
        Regression test of calwebb_spec2 pipeline performed on NIRSpec fixed-slit data
        that uses a single-slit subarray (S200B1).
        """
        input_file = self.get_data(self.test_dir,
                                   'jw84600002001_02101_00001_nrs2_rate.fits')
        step = Spec2Pipeline()
        step.save_bsub = True
        step.save_results = True
        step.resample_spec.save_results = True
        step.cube_build.save_results = True
        step.extract_1d.save_results = True
        step.run(input_file)

        outputs = [('jw84600002001_02101_00001_nrs2_cal.fits',
                    'jw84600002001_02101_00001_nrs2_cal_ref.fits'),
                   ('jw84600002001_02101_00001_nrs2_s2d.fits',
                    'jw84600002001_02101_00001_nrs2_s2d_ref.fits'),
                   ('jw84600002001_02101_00001_nrs2_x1d.fits',
                    'jw84600002001_02101_00001_nrs2_x1d_ref.fits')
                  ]
        self.compare_outputs(outputs)

@pytest.mark.bigdata
class TestNIRSpecMasterBackGround_FS(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_masterbackground', 'nrs-fs', 'truth']
    test_dir = ['test_masterbackground','nrs-fs']
    rtol = 0.0001

    def test_nirspec_masterbackground_fs_user1d(self):
        """

        Regression test of master background subtraction got NRS FS  when a user 1-D spectrum is provided.

        """
        # input file has 2-D background image added to it
        input_file = self.get_data(*self.test_dir,
                                  'nrs1_sci_bkg_cal.fits')
        # user provide 1-D background was created from the 2-D background image
        input_1dbgd_file = self.get_data(*self.test_dir,
                                         'nrs1_bkg_x1d_user.fits')

        #result = MasterBackgroundStep.call(input_file,user_background=input_1dbgd_file, 
        #save_results=True)
       
        #____________________________________________________________________________
        # Test 1 compare 4 FS 1D extracted spectra from sciene data with no background added
        # to 4 FS 1D extracted spectra from the output from MasterBackground subtraction

        # run 1-D extract on results from MasterBackground step
        #result_1d = Extract1dstep.call(result,save_results=True)
        # just to get some data to test set make up this data
        input_test= self.get_data(*self.test_dir, 
                                            'nrs1_sci_cal.fits')  
        result_1d = Extract1dStep.call(input_test,save_results=True)

        # run 1-D extract on original science data without background
        input_sci_cal_file = self.get_data(*self.test_dir, 
                                 'nrs1_sci_cal.fits')  
        # Should I instead run this off line to produce this file and make this
        # one of the input files or truth files ?
        sci_cal_1d = Extract1dStep.call(input_sci_cal_file,save_results=True)

        # Compare the FS  1D extracted data. These types of data are
        #  MultiSpec Models
        num_spec = len(result_1d.spec)
        atol = 0.01
        rtol = 0.0005
        for i in range(num_spec):
            sci_spec_1d = sci_cal_1d.spec[i].spec_table.FLUX
            result_spec_1d = result_1d.spec[i].spec_table.FLUX
            assert_allclose(sci_spec_1d,result_spec_1d,rtol=rtol,atol=atol)
        #____________________________________________________________________________
        # Test 2  compare the science MultiSlit data with no background to the
        # MultiSlit output from the masterBackground Subtraction step 
        # background subtracted science image. For FS case 
        result = datamodels.open(input_test) # override result just for testing

        input_sci_model = datamodels.open(input_sci_cal_file)
        num_slits = len(input_sci_model.slits)
        
        atol = 0.01
        rtol = 0.0005
        for i in range(num_slits):
            slit_sci = input_sci_model.slits[i].data
            slit_result = result.slits[i].data
            
            sub = slit_sci - slit_result
            sub_smo = convolve(sub,Box2DKernel(3))
            sub_smo_zero = sub_smo*0.0
            assert_allclose(sub_smo,sub_smo_zero,rtol=rtol,atol=atol)


        #____________________________________________________________________________
        # Test 3 Compare background sutracted science data (results) to truth file 
        # This data is MultiSlit data

        ref_file = self.get_data(*self.ref_loc,
                                 'nrs1_sci_cal.fits') # temp file replace when we have 
                                                      # a real truth file

        result_file = result.meta.filename
        result.save(result_file)

        
        outputs = [(result_file,ref_file)]
        self.compare_outputs(outputs)
        result.close()
