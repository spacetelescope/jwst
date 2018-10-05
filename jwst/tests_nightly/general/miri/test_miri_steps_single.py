import os
import numpy as np
from numpy.testing import utils
import pytest

from jwst import datamodels
from jwst.datamodels import ImageModel, RegionsModel, CubeModel
from jwst.stpipe import crds_client
from jwst.lib.set_telescope_pointing import add_wcs

from jwst.tests.base_test import BaseJWSTTest

from jwst.assign_wcs import AssignWcsStep
from jwst.cube_build.cube_build_step import CubeBuildStep
from jwst.linearity.linearity_step import LinearityStep
from jwst.ramp_fitting.ramp_fit_step import RampFitStep


@pytest.mark.bigdata
class TestMIRIRampFit(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_ramp_fit', 'truth']
    test_dir = 'test_ramp_fit'

    def test_ramp_fit_miri1(self):
        """
        Regression test of ramp_fit step performed on MIRI data.
        """
        input_file = self.get_data(self.test_dir, 'jw00001001001_01101_00001_MIRIMAGE_jump.fits')

        result = RampFitStep.call(input_file,
                         save_opt=True,
                         opt_name='rampfit1_opt_out.fits')
        output_file = result[0].save(path=result[0].meta.filename.replace('jump','rampfit'))
        int_output = result[1].save(path=result[1].meta.filename.replace('jump','rampfit_int'))
        result[0].close()
        result[1].close()
        
        outputs = [(output_file, 
                    'jw00001001001_01101_00001_MIRIMAGE_ramp_fit.fits'),
                   (int_output,
                    'jw00001001001_01101_00001_MIRIMAGE_int.fits'),
                   ('rampfit1_opt_out_fitopt.fits',
                    'jw00001001001_01101_00001_MIRIMAGE_opt.fits')
                  ]
        self.compare_outputs(outputs)

    def test_ramp_fit_miri2(self):
        """
        Regression test of ramp_fit step performed on MIRI data.
        """
        input_file = self.get_data(self.test_dir, 
                                   'jw80600012001_02101_00003_mirimage_jump.fits')

        result = RampFitStep.call(input_file,
                          save_opt=True,
                          opt_name='rampfit2_opt_out.fits')
             
        output_file = result[0].save(path=result[0].meta.filename.replace('jump','rampfit'))
        int_output = result[1].save(path=result[1].meta.filename.replace('jump','rampfit_int'))
        result[0].close()
        result[1].close()
        
        outputs = [(output_file, 
                    'jw80600012001_02101_00003_mirimage_ramp.fits'),
                   (int_output,
                    'jw80600012001_02101_00003_mirimage_int.fits'),
                   ('rampfit2_opt_out_fitopt.fits',
                    'jw80600012001_02101_00003_mirimage_opt.fits')
                  ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestMIRICube(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_cube_build', 'truth']
    test_dir = 'test_cube_build'

    def test_cubebuild_miri(self):
        """
        Regression test of cube_build performed on MIRI MRS data.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw10001001001_01101_00001_mirifushort_cal.fits')

        input_model = datamodels.IFUImageModel(input_file)
        result = CubeBuildStep.call(input_model, output_type='multi', save_results=True)

        paths = result.save()
        result.close()
        output_file = paths[0]

        outputs = [(output_file,
                    'jw10001001001_01101_00001_mirifushort_s3d_ref.fits',
                    ['primary','sci','err','dq','wmap']) ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestMIRILinearity(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_linearity','truth']
    test_dir ='test_linearity'

    def test_linearity_miri3(self):
        """
        Regression test of linearity step performed on MIRI data.
        """
        input_file = self.get_data(self.test_dir,
                                    'jw00001001001_01109_00001_MIRIMAGE_dark_current.fits')
        # get supplemental input
        override_file = self.get_data(self.test_dir,
                                      "lin_nan_flag_miri.fits")
        # run calibration step
        result = LinearityStep.call(input_file,
                           override_linearity=override_file)

        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        outputs = [(output_file,
                    'jw00001001001_01109_00001_MIRIMAGE_linearity.fits') ]
        self.compare_outputs(outputs)


@pytest.mark.bigdata
class TestMIRIWCSFixed(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_wcs','fixed','truth']
    test_dir = os.path.join('test_wcs','fixed')

    def test_miri_fixed_slit_wcs(self):
        """
        Regression test of creating a WCS object and doing pixel to sky transformation.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw00035001001_01101_00001_mirimage_rate.fits')
        ref_file = self.get_data(os.path.join(*self.ref_loc),
                                 'jw00035001001_01101_00001_mirimage_assign_wcs.fits')

        result = AssignWcsStep.call(input_file)
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        im = ImageModel(output_file)
        imref = ImageModel(ref_file)
        y, x = np.mgrid[:1031, :1024]
        ra, dec, lam = im.meta.wcs(x, y)
        raref, decref, lamref = imref.meta.wcs(x, y)
        utils.assert_allclose(ra, raref)
        utils.assert_allclose(dec, decref)
        utils.assert_allclose(lam, lamref)


@pytest.mark.bigdata
class TestMIRIWCSIFU(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_wcs', 'ifu', 'truth']
    test_dir = os.path.join('test_wcs', 'ifu')

    def test_miri_ifu_wcs(self):
        """
        Regression test of creating a WCS object and doing pixel to sky transformation.
        """
        input_file = self.get_data(self.test_dir,
                                   'jw00024001001_01101_00001_MIRIFUSHORT_uncal_MiriSloperPipeline.fits')
        ref_file = self.get_data(os.path.join(*self.ref_loc),
                                 'jw00024001001_01101_00001_MIRIFUSHORT_assign_wcs.fits')


        result = AssignWcsStep.call(input_file)
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        im = ImageModel(output_file)
        imref = ImageModel(ref_file)

        # Get the region file
        region = RegionsModel(crds_client.get_reference_file(im, 'regions'))

        # inputs
        shape = region.regions.shape
        y, x = np.mgrid[ : shape[0], : shape[1]]

        # Get indices where pixels == 0. These should be NaNs in the output.
        ind_zeros = region.regions == 0

        ra, dec, lam = im.meta.wcs(x, y)
        raref, decref, lamref = imref.meta.wcs(x, y)
        utils.assert_allclose(ra, raref, equal_nan=True)
        utils.assert_allclose(dec, decref, equal_nan=True)
        utils.assert_allclose(lam, lamref, equal_nan=True)

        # Test that we got NaNs at ind_zero
        assert(np.isnan(ra).nonzero()[0] == ind_zeros.nonzero()[0]).all()
        assert(np.isnan(ra).nonzero()[1] == ind_zeros.nonzero()[1]).all()

        # Test the inverse transform
        x1, y1 = im.meta.wcs.backward_transform(ra, dec, lam)
        assert(np.isnan(x1).nonzero()[0] == ind_zeros.nonzero()[0]).all()
        assert (np.isnan(x1).nonzero()[1] == ind_zeros.nonzero()[1]).all()

        # Also run a smoke test with values outside the region.
        dec[100][200] = -80
        ra[100][200] = 7
        lam[100][200] = 15

        x2, y2 = im.meta.wcs.backward_transform(ra, dec, lam)
        assert np.isnan(x2[100][200])
        assert np.isnan(x2[100][200])


@pytest.mark.bigdata
class TestMIRIWCSImage(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_wcs','image','truth']
    test_dir = os.path.join('test_wcs','image')

    def test_miri_image_wcs(self):
        """
        Regression test of creating a WCS object and doing pixel to sky transformation.
        """

        input_file = self.get_data(self.test_dir,
                                    "jw00001001001_01101_00001_MIRIMAGE_ramp_fit.fits")
        ref_file = self.get_data(*self.ref_loc,
                                 "jw00001001001_01101_00001_MIRIMAGE_assign_wcs.fits")

        result = AssignWcsStep.call(input_file)
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        im = ImageModel(output_file)
        imref = ImageModel(ref_file)
        x, y = np.mgrid[:1031, :1024]
        ra, dec = im.meta.wcs(x, y)
        raref, decref = imref.meta.wcs(x, y)
        utils.assert_allclose(ra, raref)
        utils.assert_allclose(dec, decref)


@pytest.mark.bigdata
class TestMIRIWCSSlitless(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_wcs','slitless','truth']
    test_dir = os.path.join('test_wcs','slitless')

    def test_miri_slitless_wcs(self):
        """

        Regression test of creating a WCS object and doing pixel to sky transformation.

        """
        input_file = self.get_data(self.test_dir,
                                   "jw80600012001_02101_00003_mirimage_rateints.fits")
        ref_file = self.get_data(*self.ref_loc,
                                 "jw80600012001_02101_00003_mirimage_assign_wcs.fits")

        result = AssignWcsStep.call(input_file)
        output_file = result.meta.filename
        result.save(output_file)
        result.close()

        im = CubeModel(output_file)
        imref = CubeModel(ref_file)
        x, y = np.mgrid[:1031, :1024]
        ra, dec, lam = im.meta.wcs(x, y)
        raref, decref, lamref = imref.meta.wcs(x, y)
        utils.assert_allclose(ra, raref)
        utils.assert_allclose(dec, decref)
        utils.assert_allclose(lam, lamref)


@pytest.mark.bigdata
class TESTMIRISetPointing(BaseJWSTTest):
    input_loc = 'miri'
    ref_loc = ['test_pointing', 'truth']
    test_dir = 'test_pointing'
    rtol = 0.000001

    def test_miri_setpointing(self):
        """
        Regression test of the set_telescope_pointing script on a level-1b MIRI file.
        """

        # Copy original version of file to test file, which will get overwritten by test
        orig_file = self.get_data(self.test_dir,
                                    'jw80600010001_02101_00001_mirimage_uncal_orig.fits',
                                    docopy=True  # always produce local copy
                              )
        input_file = orig_file.replace('_orig.fits','.fits')
        os.rename(orig_file, input_file)  # rename local copy

        add_wcs(input_file)

        outputs = [(input_file,
                    'jw80600010001_02101_00001_mirimage_uncal_ref.fits')]
        self.compare_outputs(outputs)
