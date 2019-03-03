import pytest
import numpy as np

from numpy.testing import assert_allclose
from gwcs.wcstools import grid_from_bounding_box
from jwst.tests.base_classes import BaseJWSTTest
from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.datamodels import ImageModel
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline
from jwst.imprint import ImprintStep
from jwst.ramp_fitting import RampFitStep
from jwst.extract_1d import Extract1dStep
from jwst.master_background import MasterBackgroundStep
from jwst.cube_build import CubeBuildStep
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
class TestNIRSpecMasterBackground_FS(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_masterbackground', 'nrs-fs', 'truth']
    test_dir = ['test_masterbackground', 'nrs-fs']
    rtol = 0.0001

    def test_nirspec_masterbackground_fs_user1d(self):
        """

        Regression test of master background subtraction for NRS FS when a user 1-D spectrum is provided.

        """
        # input file has 2-D background image added to it
        input_file = self.get_data(*self.test_dir,
                                    'nrs_sci+bkg_cal.fits')
        # user provided 1-D background was created from the 2-D background image
        input_1dbkg_file = self.get_data(*self.test_dir,
                                          'nrs_bkg_user_clean_x1d.fits')

        result = MasterBackgroundStep.call(input_file,
                                           user_background=input_1dbkg_file)
        # _________________________________________________________________________
        # Test 1 compare 4 FS 1D extracted spectra from science data with
        # no background added to 4 FS 1D extracted spectra from the output
        # from MasterBackground subtraction

        # run 1-D extract on results from MasterBackground step
        result_1d = Extract1dStep.call(result)

        # run 1-D extract on original science data without background
        input_sci_cal_file = self.get_data(*self.test_dir,
                                            'nrs_sci_cal.fits')
        # run 1D extraction on original science data without background
        sci_cal_1d = Extract1dStep.call(input_sci_cal_file)

        # Compare the FS  1D extracted data. These types of data are
        #  MultiSpec Models.
        input_sci = datamodels.open(input_sci_cal_file)
        num_spec = len(sci_cal_1d.spec)

        # the user 1D spectrum may not cover the entire wavelength range of the
        # science data.  Find the wavelength range of user 1-D spectra
        input_1dbkg_1d = datamodels.open(input_1dbkg_file)
        user_wave = input_1dbkg_1d.spec[0].spec_table['wavelength']
        user_flux = input_1dbkg_1d.spec[0].spec_table['flux']
        user_wave_valid = np.where(user_flux > 0)
        min_user_wave = np.amin(user_wave[user_wave_valid])
        max_user_wave = np.amax(user_wave[user_wave_valid])
        input_1dbkg_1d.close()

        for i in range(num_spec):
            # ______________________________________________________________________
            # Test 1 compare extracted spectra data from the science data
            # to extracted spectra from the output
            # from MasterBackground subtraction.
            sci_spec_1d = sci_cal_1d.spec[i].spec_table['flux']
            sci_wave = sci_cal_1d.spec[i].spec_table['wavelength']
            result_spec_1d = result_1d.spec[i].spec_table['flux']

            # find the waverange covered by both user 1-D and science slit
            sci_wave_valid = np.where(sci_spec_1d > 0)
            min_wave = np.amin(sci_wave[sci_wave_valid])
            max_wave = np.amax(sci_wave[sci_wave_valid])
            if min_user_wave > min_wave:
                min_wave = min_user_wave
            if max_user_wave < max_wave:
                max_wave = max_user_wave

            sub_spec = sci_spec_1d - result_spec_1d
            valid = np.where(np.logical_and(sci_wave > min_wave, sci_wave < max_wave))
            sub_spec = sub_spec[valid]
            mean_sub = np.absolute(np.nanmean(sub_spec))
            atol = 4.0
            rtol = 0.02
            assert_allclose(mean_sub, 0, atol=atol)
            # ______________________________________________________________________
            # Test 2  compare the science  data with no background
            # to the output from the masterBackground Subtraction step
            # background subtracted science image.
            bb = input_sci.slits[i].meta.wcs.bounding_box
            x, y = grid_from_bounding_box(bb)
            ra, dec, lam = input_sci.slits[i].meta.wcs(x, y)
            valid = np.isfinite(lam)

            sci = input_sci.slits[i].data
            result_slit = result.slits[i].data

            # check for outliers in the science image that could cause test
            # to fail. These could be cosmic ray hits or other yeck that
            # messes up the science data - ignores these cases
            sci_mean = np.nanmean(sci[valid])
            sci_std = np.nanstd(sci[valid])
            upper = sci_mean + sci_std*5.0
            lower = sci_mean - sci_std*5.0
            mask_clean = np.logical_and(sci[valid] < upper, sci[valid] > lower)

            # for this slit subtract the background subtracted data from
            # the science data with no background added
            sub = result_slit - sci
            # do not look at outliers - they confuse things
            sub_valid = sub[valid]
            mean_sub = np.mean(sub_valid[mask_clean])
            atol = 0.5
            assert_allclose(np.absolute(mean_sub), 0, atol=atol)
        input_sci.close()
        # ______________________________________________________________________
        # Test 3 Compare background sutracted science data (results)
        #  to a truth file. This data is MultiSlit data
        result_file = result.meta.filename
        result.save(result_file)
        result.close()
        ref_file = self.get_data(*self.ref_loc,
                                  'nrs_sci+bkg_masterbackgroundstep.fits')

        outputs = [(result_file, ref_file)]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestNIRSpecMasterBackground_IFU(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_masterbackground', 'nrs-ifu', 'truth']
    test_dir = ['test_masterbackground', 'nrs-ifu']
    rtol = 0.0001

    def test_nirspec_masterbackground_ifu_user1d(self):
        """

        Regression test of master background subtraction for NRS IFU when a user 1-D spectrum is provided.

        """
        # input file has 2-D background image added to it

        input_file = self.get_data(*self.test_dir,
                                    'prism_sci_bkg_cal.fits')
        # user provide 1-D background was created from the 2-D background image
        input_1dbkg_file = self.get_data(*self.test_dir,
                                          'prism_bkg_x1d.fits')

        result = MasterBackgroundStep.call(input_file,
                                           user_background=input_1dbkg_file)

        # _________________________________________________________________________
        # Test 1 compare extracted spectra data with
        # no background added to extracted spectra from the output
        # from MasterBackground subtraction. First cube_build has to be run
        # on the data.

        result_s3d = CubeBuildStep.call(result)
        # run 1-D extract on results from MasterBackground step
        result_1d = Extract1dStep.call(result_s3d, subtract_background=False)

        # run 1-D extract on original science data without background
        input_sci_cal_file = self.get_data(*self.test_dir,
                                            'prism_sci_cal.fits')
        # find the 1D extraction of this file
        # find Extract1dStep on sci data to use the same version of
        # this rountine run on both files  rather than running it off line
        # and having the 1-D extracted science file stored as input to step

        sci_s3d = CubeBuildStep.call(input_sci_cal_file)
        sci_1d = Extract1dStep.call(sci_s3d, subtract_background=False)

        # read in the valid wavelengths of the user-1d
        input_1d_bkg_model = datamodels.open(input_1dbkg_file)
        user_wave = input_1d_bkg_model.spec[0].spec_table['wavelength']
        user_flux = input_1d_bkg_model.spec[0].spec_table['net']
        user_wave_valid = np.where(user_flux > 0)
        min_user_wave = np.amin(user_wave[user_wave_valid])
        max_user_wave = np.amax(user_wave[user_wave_valid])
        input_1d_bkg_model.close()
        # find the waverange covered by both user and science
        sci_spec_1d = sci_1d.spec[0].spec_table['net']
        sci_spec_wave = sci_1d.spec[0].spec_table['wavelength']

        result_spec_1d = result_1d.spec[0].spec_table['net']

        sci_wave_valid = np.where(sci_spec_1d > 0)
        min_wave = np.amin(sci_spec_wave[sci_wave_valid])
        max_wave = np.amax(sci_spec_wave[sci_wave_valid])
        if min_user_wave > min_wave:
            min_wave = min_user_wave
        if max_user_wave < max_wave:
            max_wave = max_user_wave

        sub_spec = sci_spec_1d - result_spec_1d
        valid = np.where(np.logical_and(sci_spec_wave > min_wave, sci_spec_wave < max_wave))
        sub_spec = sub_spec[valid]
        sub_spec = sub_spec[1:-2] #  endpoints are wacky

        mean_sub = np.absolute(np.nanmean(sub_spec))
        atol = 5.0
        assert_allclose(mean_sub, 0, atol=atol)
        # ______________________________________________________________________
        # Test 2  compare the science  data with no background
        # to the output from the masterBackground Subtraction step
        # background subtracted science image.
        input_sci_model = datamodels.open(input_sci_cal_file)

        # We don't want the slices gaps to impact the statisitic
        # loop over the 30 Slices
        for i in range(30):
            slice_wcs = nirspec.nrs_wcs_set_input(input_sci_model, i)
            x, y = grid_from_bounding_box(slice_wcs.bounding_box)
            ra, dec, lam = slice_wcs(x, y)
            valid = np.isfinite(lam)
            result_slice_region = result.data[y.astype(int), x.astype(int)]
            sci_slice_region = input_sci_model.data[y.astype(int),
                                                    x.astype(int)]
            sci_slice = sci_slice_region[valid]
            result_slice = result_slice_region[valid]
            sub = result_slice - sci_slice

            # check for outliers in the science image
            sci_mean = np.nanmean(sci_slice)
            sci_std = np.nanstd(sci_slice)
            upper = sci_mean + sci_std*5.0
            lower = sci_mean - sci_std*5.0
            mask_clean = np.logical_and(sci_slice < upper, sci_slice > lower)

            sub_mean = np.absolute(np.nanmean(sub[mask_clean]))
            atol = 2.0
            assert_allclose(sub_mean, 0, atol=atol)
        # ______________________________________________________________________
        # Test 3 Compare background sutracted science data (results)
        #  to a truth file. This data is MultiSlit data

        input_sci_model.close()
        result_file = result.meta.filename
        result.save(result_file)
        result.close()
        ref_file = self.get_data(*self.ref_loc,
                                  'prism_sci_bkg_masterbackgroundstep.fits')

        outputs = [(result_file, ref_file)]
        self.compare_outputs(outputs)
        input_sci_model.close()


@pytest.mark.bigdata
class TestNIRSpecMasterBackground_MOS(BaseJWSTTest):
    input_loc = 'nirspec'
    ref_loc = ['test_masterbackground', 'nrs-mos', 'truth']
    test_dir = ['test_masterbackground', 'nrs-mos']
    rtol = 0.0001

    def test_nirspec_masterbackground_mos_user1d(self):
        """

        Regression test of master background subtraction for NRS MOS when a user 1-D spectrum is provided.

        """
        # input file has 2-D background image added to it
        input_file = self.get_data(*self.test_dir,
                                    'nrs_mos_sci+bkg_cal.fits')
        # user provide 1-D background was created from the 2-D background image
        input_1dbkg_file = self.get_data(*self.test_dir,
                                          'nrs_mos_bkg_x1d.fits')

        result = MasterBackgroundStep.call(input_file,
                                           user_background=input_1dbkg_file)
        # _________________________________________________________________________
        # One of out tests is to compare the 1-D extracted spectra from
        # the science image (no background added) and the masterbackground subtracted
        # data. Run extract1d on both of these files

        # run 1-D extract on results from MasterBackground step
        result_1d = Extract1dStep.call(result)

        # run 1-D extract on original science data without background
        input_sci_cal_file = self.get_data(*self.test_dir,
                                            'nrs_mos_sci_cal.fits')

        input_sci = datamodels.open(input_sci_cal_file)
        sci_cal_1d = Extract1dStep.call(input_sci)
        num_spec = len(sci_cal_1d.spec)

        # the user 1D spectrum may not cover the entire wavelength range of the
        # science data.  Find the wavelength range of user 1-D spectra

        input_1dbkg_1d = datamodels.open(input_1dbkg_file)
        user_wave = input_1dbkg_1d.spec[0].spec_table['wavelength']
        user_flux = input_1dbkg_1d.spec[0].spec_table['flux']
        user_wave_valid = np.where(user_flux > 0)
        min_user_wave = np.amin(user_wave[user_wave_valid])
        max_user_wave = np.amax(user_wave[user_wave_valid])
        input_1dbkg_1d.close()

        # loop over each slit and perform 2 tests on each slit
        for i in range(num_spec):
            # ______________________________________________________________________
            # Test 1 compare extracted spectra data from the science data
            # to extracted spectra from the output
            # from MasterBackground subtraction.
            sci_spec_1d = sci_cal_1d.spec[i].spec_table['flux']
            sci_wave = sci_cal_1d.spec[i].spec_table['wavelength']
            result_spec_1d = result_1d.spec[i].spec_table['flux']

            # find the waverange covered by both user 1-D and science slit
            sci_wave_valid = np.where(sci_spec_1d > 0)
            min_wave = np.amin(sci_wave[sci_wave_valid])
            max_wave = np.amax(sci_wave[sci_wave_valid])
            if min_user_wave > min_wave:
                min_wave = min_user_wave
            if max_user_wave < max_wave:
                max_wave = max_user_wave

            sub_spec = sci_spec_1d - result_spec_1d
            valid = np.where(np.logical_and(sci_wave > min_wave, sci_wave < max_wave))
            sub_spec = sub_spec[valid]
            mean_sub = np.nanmean(sub_spec)
            atol = 3.0
            assert_allclose(mean_sub, 0, atol=atol)
            # ______________________________________________________________________
            # Test 2  compare the science  data with no background
            # to the output from the masterBackground Subtraction step
            # background subtracted science image.
            bb = input_sci.slits[i].meta.wcs.bounding_box
            x, y = grid_from_bounding_box(bb)
            ra, dec, lam = input_sci.slits[i].meta.wcs(x, y)
            valid = np.isfinite(lam)

            sci = input_sci.slits[i].data
            result_slit = result.slits[i].data

            # check for outliers in the science image that could cause test
            # to fail. These could be cosmic ray hits or other yeck that
            # messes up the science data - ignores these cases
            sci_mean = np.nanmean(sci[valid])
            sci_std = np.nanstd(sci[valid])
            upper = sci_mean + sci_std*5.0
            lower = sci_mean - sci_std*5.0
            mask_clean = np.logical_and(sci[valid] < upper, sci[valid] > lower)

            # for this slit subtract the background subtracted data from
            # the science data with no background added
            sub = result_slit - sci
            # do not look at outliers - they confuse things
            sub_valid = sub[valid]
            mean_sub = np.mean(sub_valid[mask_clean])
            atol = 0.1
            assert_allclose(np.absolute(mean_sub), 0, atol=atol)
        # ______________________________________________________________________
        # Test 3 Compare background sutracted science data (results)
        #  to a truth file. This data is MultiSlit data

        result_file = result.meta.filename
        result.save(result_file)
        result.close()
        ref_file = self.get_data(*self.ref_loc,
                                  'nrs_mos_sci+bkg_masterbackgroundstep.fits')

        outputs = [(result_file, ref_file)]
        self.compare_outputs(outputs)
        input_sci.close()
