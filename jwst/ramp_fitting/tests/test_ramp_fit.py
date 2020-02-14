import pytest
import numpy as np

from jwst.ramp_fitting.ramp_fit import ramp_fit
from jwst.datamodels import dqflags
from jwst.datamodels import RampModel
from jwst.datamodels import GainModel, ReadnoiseModel

# single group intergrations fail in the GLS fitting
# so, keep the two method test separate and mark GLS test as
# expected to fail.  Needs fixing, but the fix is not clear
# to me. [KDG - 19 Dec 2018]

def test_one_group_small_buffer_fit_ols():
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=1,gain=1,readnoise=10)
    model1.data[0, 0, 50, 50] = 10.0
    slopes = ramp_fit(model1, 512, True, rnModel, gain, 'OLS', 'optimal')
    np.testing.assert_allclose(slopes[0].data[50, 50],10.0, 1e-6)

@pytest.mark.xfail(reason="GLS code does not [yet] handle single group integrations.")
def test_one_group_small_buffer_fit_gls():
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=1,gain=1,readnoise=10)
    model1.data[0, 0, 50, 50] = 10.0
    slopes = ramp_fit(model1, 512, True, rnModel, gain, 'GLS', 'optimal')
    np.testing.assert_allclose(slopes[0].data[50, 50],10.0, 1e-6)

def test_one_group_two_ints_fit_ols():
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=1,gain=1,readnoise=10,nints=2)
    model1.data[0, 0, 50, 50] = 10.0
    model1.data[1, 0, 50, 50] = 12.0
    slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
    np.testing.assert_allclose(slopes[0].data[50, 50],11.0, 1e-6)

@pytest.mark.xfail(reason="GLS code does not [yet] handle single group integrations.")
def test_one_group_two_ints_fit_gls():
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=1,gain=1,readnoise=10,nints=2)
    model1.data[0, 0, 50, 50] = 10.0
    model1.data[1, 0, 50, 50] = 12.0
    slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'GLS', 'optimal')
    np.testing.assert_allclose(slopes[0].data[50, 50],11.0, 1e-6)

# tests that apply to both 'ols' and 'gls' are in the TestMethods class so
# that both can use the parameterized 'method'

@pytest.mark.parametrize("method", ['OLS', 'GLS'])
class TestMethods:

    def test_nocrs_noflux(self, method):
        # all pixel values are zero. So slope should be zero
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
        slopes = ramp_fit(model1, 60000, False, rnModel, gain, method, 'optimal')
        assert(0 == np.max(slopes[0].data))
        assert(0 == np.min(slopes[0].data))

    def test_nocrs_noflux_firstrows_are_nan(self, method):
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
        model1.data[0, : , 0:12, :] = np.nan
        slopes = ramp_fit(model1, 60000, False, rnModel, gain, 'OLS', 'optimal')
        assert(0 == np.max(slopes[0].data))
        assert(0 == np.min(slopes[0].data))

    @pytest.mark.xfail(reason="Fails, without frame_time it doesn't work")
    def test_error_when_frame_time_not_set(self, method):
            #all pixel values are zero. So slope should be zero
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
        model1.meta.exposure.frame_time = None
        slopes = ramp_fit(model1, 64000, False, rnModel, gain, 'OLS', 'optimal')
        assert(0 == np.max(slopes[0].data))
        assert(0 == np.min(slopes[0].data))

    def test_five_groups_two_ints_Poisson_noise_only(self, method):
        grouptime=3.0
        ingain = 2000
        inreadnoise =7
        ngroups=5
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                              gain=ingain, readnoise=inreadnoise,deltatime=grouptime, nints=2)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        model1.data[1, 0, 50, 50] = 10.0
        model1.data[1, 1, 50, 50] = 15.0
        model1.data[1, 2, 50, 50] = 25.0
        model1.data[1, 3, 50, 50] = 33.0
        model1.data[1, 4, 50, 50] = 160.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        out_slope=slopes[0].data[50, 50]
        deltaDN1 = 50
        deltaDN2 = 150
        np.testing.assert_allclose(out_slope, (deltaDN1 + deltaDN2)/2.0, 75.0, 1e-6)

    def test_ngroups_doesnot_match_cube_size(self, method):
            #all pixel values are zero. So slope should be zero
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
        model1.meta.exposure.ngroups = 11
        slopes = ramp_fit(model1, 64000, False, rnModel, gain, 'OLS', 'optimal')
        assert(0 == np.max(slopes[0].data))
        assert(0 == np.min(slopes[0].data))

    def test_bad_gain_values(self, method):
            #all pixel values are zero. So slope should be zero
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
        model1.meta.exposure.ngroups = 11
        gain.data[10,10] = -10
        gain.data[20,20] = np.nan
        slopes = ramp_fit(model1, 64000, False, rnModel, gain, 'OLS', 'optimal')
        assert(0 == np.max(slopes[0].data))
        assert(0 == np.min(slopes[0].data))
        assert slopes[0].dq[10, 10] == 524288 + 1
        assert slopes[0].dq[20, 20] == 524288 + 1

    def test_subarray_5groups(self, method):
            #all pixel values are zero. So slope should be zero
        model1, gdq, rnModel, pixdq, err, gain = setup_subarray_inputs(ngroups=5,
            subxstart=10, subystart=20, subxsize=5, subysize=15, readnoise=50)
        model1.meta.exposure.ngroups = 11
        model1.data[0, 0, 12, 1] = 10.0
        model1.data[0, 1, 12, 1] = 15.0
        model1.data[0, 2, 12, 1] = 25.0
        model1.data[0, 3, 12, 1] = 33.0
        model1.data[0, 4, 12, 1] = 60.0
        xvalues = np.arange(5)*1.0
        yvalues = np.array([10,15,25,33,60])
        coeff = np.polyfit(xvalues, yvalues, 1)
        slopes = ramp_fit(model1, 64000, False, rnModel, gain, 'OLS', 'optimal')
        np.testing.assert_allclose(slopes[0].data[12,1],coeff[0],1e-6)

    def test_simple_ramp(self, method):
        #Here given a 10 group ramp with an exact slope of 20/group. The output slope should be 20.
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=10, deltatime=3)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 30.0
        model1.data[0, 2, 50, 50] = 50.0
        model1.data[0, 3, 50, 50] = 70.0
        model1.data[0, 4, 50, 50] = 90.0
        model1.data[0, 5, 50, 50] = 110.0
        model1.data[0, 6, 50, 50] = 130.0
        model1.data[0, 7, 50, 50] = 150.0
        model1.data[0, 8, 50, 50] = 170.0
        model1.data[0, 9, 50, 50] = 190.0
        slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
        # take the ratio of the slopes to get the relative error
        np.testing.assert_allclose(slopes[0].data[50, 50], (20.0/3), 1e-6)

    def test_read_noise_only_fit(self, method):
            #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,readnoise=50)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        xvalues = np.arange(5)*1.0
        yvalues = np.array([10,15,25,33,60])
        coeff = np.polyfit(xvalues, yvalues, 1)
        np.testing.assert_allclose(slopes[0].data[50, 50], coeff[0], 1e-6)

    def test_photon_noise_only_fit(self, method):
            #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,gain=1000,readnoise=1)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,4,50,50] - model1.data[0,0,50,50])/ 4.0
        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-2)

    def test_photon_noise_only_bad_last_frame(self, method):
            #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,gain=1000,readnoise=1)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        model1.groupdq[0,4,:,:] = dqflags.group['DO_NOT_USE']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,3,50,50] - model1.data[0,0,50,50])/ 3.0
        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-2)

    @pytest.mark.xfail(reason="Fails, bad last frame yields only one good one. \
        This should not every happen. When ngroups==2 the last frame doesn't get flagged.")
    def test_photon_noise_only_bad_last_frame_two_groups(self, method):
            #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=2,gain=1000,readnoise=1)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.groupdq[0,1,:,:] = dqflags.group['DO_NOT_USE']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,1,50,50] - model1.data[0,0,50,50])/ 1.0
        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-6)

    def test_photon_noise_with_unweighted_fit(self, method):
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,gain=1000,readnoise=1)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'unweighted')
        # cds_slope = (model1.data[0,4,500,500] - model1.data[0,0,500,500])/ 4.0
        xvalues = np.arange(5)*1.0
        yvalues = np.array([10,15,25,33,60])
        coeff = np.polyfit(xvalues, yvalues, 1)
        np.testing.assert_allclose(slopes[0].data[50, 50], coeff[0], 1e-6)

    def test_two_groups_fit(self, method):
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=2,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 0, 50, 51] = 20.0
        model1.data[0, 1, 50, 51] = 60.0
        model1.data[0, 0, 50, 52] = 200.0
        model1.data[0, 1, 50, 52] = 600.0
        model1.meta.exposure.drop_frames1 = 0
        #2nd group is saturated
        model1.groupdq[0,1,50,51]=dqflags.group['SATURATED']
        #1st group is saturated
        model1.groupdq[0,0,50,52]=dqflags.group['SATURATED']
        model1.groupdq[0,1,50,52]=dqflags.group['SATURATED'] #should not be set this way
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,1,50,50] - model1.data[0,0,50,50])
        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-6)
        #expect SATURATED
        assert slopes[0].dq[50, 51] == 2 # is there a better way to do this test?
        #expect SATURATED since 1st group is Saturated
        assert slopes[0].dq[50, 52] == 2 # is there a better way to do this test?

    def test_four_groups_oneCR_orphangroupatend_fit(self, method):
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=4,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 20.0
        model1.data[0, 3, 50, 50] = 145.0
        model1.groupdq[0,3,50,50]=dqflags.group['JUMP_DET']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,1,50,50] - model1.data[0,0,50,50])

        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-6)

    #@pytest.mark.skip(reason="not using now")
    def test_four_groups_two_CRs_at_end(self, method):
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=4,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 145.0
        model1.groupdq[0,2,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,3,50,50]=dqflags.group['JUMP_DET']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        cds_slope = (model1.data[0,1,50,50] - model1.data[0,0,50,50])
        np.testing.assert_allclose(slopes[0].data[50, 50], cds_slope, 1e-6)

    def test_four_groups_four_CRs(self, method):
        #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=4,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 145.0
        model1.groupdq[0,0,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,1,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,2,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,3,50,50]=dqflags.group['JUMP_DET']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        np.testing.assert_allclose(slopes[0].data[50, 50], 0,1e-6)

    def test_four_groups_three_CRs_at_end(self, method):
        #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=4,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 145.0
        model1.groupdq[0,1,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,2,50,50]=dqflags.group['JUMP_DET']
        model1.groupdq[0,3,50,50]=dqflags.group['JUMP_DET']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        expected_slope=10.0
        np.testing.assert_allclose(slopes[0].data[50, 50],expected_slope, 1e-6)

    def test_four_groups_CR_causes_orphan_1st_group(self, method):
        #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=4,gain=.01,readnoise=10000)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 125.0
        model1.data[0, 2, 50, 50] = 145.0
        model1.data[0, 3, 50, 50] = 165.0
        model1.groupdq[0,1,50,50]=dqflags.group['JUMP_DET']
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        expected_slope=20.0
        np.testing.assert_allclose(slopes[0].data[50, 50],expected_slope, 1e-6)

    def test_one_group_fit(self, method):
        #
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=1,gain=1,readnoise=10)
        model1.data[0, 0, 50, 50] = 10.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        np.testing.assert_allclose(slopes[0].data[50, 50],10.0, 1e-6)

    def test_two_groups_unc(self, method):
        grouptime=3.0
        deltaDN = 5
        ingain = 2
        inreadnoise =10
        ngroups=2
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
            gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 10.0 + deltaDN
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        # delta_electrons = deltaDN * ingain
        single_sample_readnoise = inreadnoise/np.sqrt(2)
        np.testing.assert_allclose(slopes[0].var_poisson[50,50],
            ((deltaDN/ingain)/grouptime**2), 1e-6)
        np.testing.assert_allclose(slopes[0].var_rnoise[50,50],
            (inreadnoise**2/grouptime**2), 1e-6)
        np.testing.assert_allclose(slopes[0].var_rnoise[50,50],
            (12*single_sample_readnoise**2/(ngroups*(ngroups**2 - 1)*grouptime**2)), 1e-6)
        np.testing.assert_allclose(slopes[0].err[50,50],
            (np.sqrt((deltaDN/ingain)/grouptime**2+(inreadnoise**2/grouptime**2))), 1e-6)

    def test_five_groups_unc(self, method):
        grouptime=3.0
        # deltaDN = 5
        ingain = 2
        inreadnoise =7
        ngroups=5
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
            gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
        model1.data[0, 0, 50, 50] = 10.0
        model1.data[0, 1, 50, 50] = 15.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 33.0
        model1.data[0, 4, 50, 50] = 60.0
        slopes = ramp_fit(model1, 1024*30000., True, rnModel, gain, 'OLS', 'optimal')
        # out_slope=slopes[0].data[500, 500]
        median_slope=np.median(np.diff(model1.data[0,:,50,50]))/grouptime
        # deltaDN = 50
        delta_time = (ngroups - 1) * grouptime
        # delta_electrons = median_slope * ingain *delta_time
        single_sample_readnoise = np.float64(inreadnoise/np.sqrt(2))
        np.testing.assert_allclose(slopes[0].var_poisson[50,50],
            ((median_slope)/(ingain*delta_time)), 1e-6)
        np.testing.assert_allclose(slopes[0].var_rnoise[50,50],
            (12 * single_sample_readnoise**2/(ngroups * (ngroups**2 - 1) * grouptime**2)),  1e-6)
        np.testing.assert_allclose(slopes[0].err[50,50],
            np.sqrt(slopes[0].var_poisson[50,50]  + slopes[0].var_rnoise[50,50] ),  1e-6)

    def test_oneCR_10_groups_combination(self, method):
        grouptime=3.0
        # deltaDN = 5
        ingain = 200 # use large gain to show that Poisson noise doesn't affect the recombination
        inreadnoise = np.float64(7)
        ngroups=10
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
            gain=ingain, readnoise=inreadnoise, deltatime=grouptime)
        # two segments perfect fit, second segment has twice the slope
        model1.data[0, 0, 50, 50] = 15.0
        model1.data[0, 1, 50, 50] = 20.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 30.0
        model1.data[0, 4, 50, 50] = 35.0
        model1.data[0, 5, 50, 50] = 140.0
        model1.data[0, 6, 50, 50] = 150.0
        model1.data[0, 7, 50, 50] = 160.0
        model1.data[0, 8, 50, 50] = 170.0
        model1.data[0, 9, 50, 50] = 180.0
        model1.groupdq[0,5,50,50]=dqflags.group['JUMP_DET']
        slopes, int_model, opt_model, gls_opt_model = ramp_fit(model1,
            1024*30000.,  True, rnModel, gain, 'OLS', 'optimal')
        segment_groups  = 5
        single_sample_readnoise = np.float64(inreadnoise / np.sqrt(2))
        #check that the segment variance is as expected
        np.testing.assert_allclose(opt_model.var_rnoise[0,0,50,50],
            (12.0 * single_sample_readnoise**2/(segment_groups * (segment_groups**2 - 1) * grouptime**2)), rtol=1e-6)
        # check the combined slope is the average of the two segments since they have the same number of groups
        np.testing.assert_allclose(slopes.data[50, 50], 2.5,rtol=1e-5)
        #check that the slopes of the two segments are correct
        np.testing.assert_allclose(opt_model.slope[0,0,50, 50], 5/3.0,rtol=1e-5)
        np.testing.assert_allclose(opt_model.slope[0,1,50, 50], 10/3.0,rtol=1e-5)

    def test_oneCR_10_groups_combination_noisy2ndSegment(self, method):
        grouptime=3.0
        # deltaDN = 5
        ingain = 200 # use large gain to show that Poisson noise doesn't affect the recombination
        inreadnoise =7
        ngroups=10
        model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                              gain=ingain, readnoise=inreadnoise,deltatime=grouptime)
        # two segments perfect fit, second segment has twice the slope
        model1.data[0, 0, 50, 50] = 15.0
        model1.data[0, 1, 50, 50] = 20.0
        model1.data[0, 2, 50, 50] = 25.0
        model1.data[0, 3, 50, 50] = 30.0
        model1.data[0, 4, 50, 50] = 35.0
        model1.data[0, 5, 50, 50] = 135.0
        model1.data[0, 6, 50, 50] = 155.0
        model1.data[0, 7, 50, 50] = 160.0
        model1.data[0, 8, 50, 50] = 168.0
        model1.data[0, 9, 50, 50] = 180.0
        model1.groupdq[0,5,50,50]=dqflags.group['JUMP_DET']
        slopes, int_model, opt_model, gls_opt_model= ramp_fit(model1, 1024*30000.,  True, rnModel, gain, 'OLS', 'optimal')
        avg_slope = (opt_model.slope[0,0,50,50] + opt_model.slope[0,1,50,50])/2.0
            #even with noiser second segment, final slope should be just the average since they have the same number of groups
        np.testing.assert_allclose(slopes.data[50, 50], avg_slope,rtol=1e-5)


def test_twenty_groups_two_segments():
    '''Test to verify weighting of multiple segments in combination:
        a) gdq all 0 ; b) 1 CR (2 segments) c) 1 CR then SAT (2 segments)
    '''
    (ngroups, nints, nrows, ncols, deltatime) = (20, 1, 1, 3, 6.)
    model1, gdq, rnModel, pixdq, err, gain = setup_small_cube( ngroups,
        nints, nrows, ncols, deltatime)

    # a) ramp having gdq all 0
    model1.data[0,:,0,0] = np.arange(ngroups)*10. + 30.

    # b) ramp having 1 CR at group 15; 2 segments
    model1.data[0,:,0,1] = np.arange(ngroups)*10. + 50.
    gdq[0,15,0,1] = dqflags.group['JUMP_DET']
    model1.data[0,15:,0,1] += 1000.

    # c) ramp having 1 CR at group 2; SAT starting in group 15
    model1.data[0,:,0,2] = np.arange(ngroups)*10. + 70.
    gdq[0,2,0,2] =  dqflags.group['JUMP_DET']
    model1.data[0,2:,0,2] += 2000.
    gdq[0,15:,0,2] = dqflags.group['SATURATED']
    model1.data[0,15:,0,2] = 25000.

    new_mod, int_model, opt_model, gls_opt_model = ramp_fit(model1, 1024*30000.,
        True, rnModel, gain, 'OLS', 'optimal')

    # Check some PRI & OPT output arrays
    np.testing.assert_allclose( new_mod.data, 10./deltatime, rtol=1E-4 )

    wh_data = opt_model.slope != 0. # only test existing segments
    np.testing.assert_allclose(opt_model.slope[wh_data], 10./deltatime, rtol=1E-4)
    np.testing.assert_allclose(opt_model.yint[0,0,0,:], model1.data[0,0,0,:], rtol=1E-5)
    np.testing.assert_allclose(opt_model.pedestal[0,0,:],
        model1.data[0,0,0,:]-10., rtol=1E-5)


def test_miri_all_sat():
    ''' Test of all groups in all integrations being saturated; all output arrays
        (particularly variances) should be 0.
    '''
    (ngroups, nints, nrows, ncols, deltatime) = (3, 2, 2, 2, 6.)
    model1, gdq, rnModel, pixdq, err, gain = setup_small_cube(ngroups,
        nints, nrows, ncols, deltatime)

    model1.groupdq[:, :, :, :] = dqflags.group['SATURATED']

    new_mod, int_model, opt_model, gls_opt_model = ramp_fit(model1, 1024*30000.,
        True, rnModel, gain, 'OLS', 'optimal')

    # Check PRI output arrays
    np.testing.assert_allclose( new_mod.data, 0.0, atol=1E-6  )
    np.testing.assert_allclose( new_mod.err, 0.0, atol=1E-6  )
    np.testing.assert_allclose( new_mod.var_poisson, 0.0, atol=1E-6  )
    np.testing.assert_allclose( new_mod.var_rnoise, 0.0, atol=1E-6  )

    # Check INT output arrays
    np.testing.assert_allclose( int_model.data, 0.0, atol=1E-6  )
    np.testing.assert_allclose( int_model.err, 0.0, atol=1E-6  )
    np.testing.assert_allclose( int_model.var_poisson, 0.0, atol=1E-6  )
    np.testing.assert_allclose( int_model.var_rnoise, 0.0, atol=1E-6  )

    # Check OPT output arrays
    np.testing.assert_allclose( opt_model.slope, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.var_poisson, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.var_rnoise, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.sigslope, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.yint, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.sigyint, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.pedestal, 0.0, atol=1E-6  )
    np.testing.assert_allclose( opt_model.weights, 0.0, atol=1E-6  )


def test_miri_first_last():
    ''' This is a test of whether ramp fitting correctly handles having all 0th
        group dq flagged as DO_NOT_USE, and all final group dq flagged as
        DO_NOT_USE for MIRI data.  For 1 pixel ([1,1]) the 1st (again, 0-based)
        group is flagged as a CR.  For such a ramp, the code removes the CR-flag
        from such a CR-affected 1st group; so if it initially was 4 it is reset
        to 0 ("good"), in which case it's as if that CR was not even there.
    '''
    (ngroups, nints, nrows, ncols, deltatime) = (10, 1, 2, 2, 3.)
    model1, gdq, rnModel, pixdq, err, gain = setup_small_cube(ngroups,
        nints, nrows, ncols, deltatime)

    # Make smooth ramps having outlier SCI values in the 0th and final groups
    #   to reveal if they are included in the fit (they shouldn't be, as those
    #   groups are flagged as DO_NOT_USE)
    model1.data[0,:,0,0] = np.array([-200., 30., 40., 50., 60., 70., 80.,
        90., 100., -500.], dtype=np.float32)
    model1.data[0,:,0,1] = np.array([-300.,  80., 90., 100., 110., 120.,
        130., 140., 150., -600.], dtype=np.float32)
    model1.data[0,:,1,0] = np.array([-200., 40., 50., 60., 70., 80., 90.,
        100., 110., 900.], dtype=np.float32)
    model1.data[0,:,1,1] = np.array([-600., 140., 150., 160., 170., 180.,
        190., 200., 210., -400.], dtype=np.float32)

    # For all pixels, set gdq for 0th and final groups to DO_NOT_USE
    model1.groupdq[:,0,:,:] = dqflags.group['DO_NOT_USE']
    model1.groupdq[:,-1,:,:] = dqflags.group['DO_NOT_USE']

    # Put CR in 1st (0-based) group
    model1.groupdq[0,1,1,1] = dqflags.group['JUMP_DET']

    new_mod, int_model, opt_model, gls_opt_model = ramp_fit(model1, 1024*30000.,
        True, rnModel, gain, 'OLS', 'optimal')

    np.testing.assert_allclose( new_mod.data, 10./3., rtol=1E-5)


def setup_small_cube(ngroups=10, nints=1, nrows=2, ncols=2, deltatime=10.,
        gain=1., readnoise =10.):
    '''Create input MIRI datacube having the specified dimensions
    '''
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    pixdq = np.zeros(shape=(nrows, ncols), dtype=np.int32)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.uint8)
    model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq)

    model1.meta.instrument.name='MIRI'
    model1.meta.instrument.detector='MIRIMAGE'
    model1.meta.instrument.filter='F480M'
    model1.meta.observation.date='2015-10-13'
    model1.meta.exposure.type='MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name='FULL'

    model1.meta.subarray.xstart = 1
    model1.meta.subarray.ystart = 1
    model1.meta.subarray.xsize = ncols
    model1.meta.subarray.ysize = nrows

    model1.meta.exposure.frame_time =deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name='MIRI'

    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = ncols
    gain.meta.subarray.ysize = nrows

    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name='MIRI'

    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = ncols
    rnModel.meta.subarray.ysize = nrows

    return model1, gdq, rnModel, pixdq, err, gain


#Need test for multi-ints near zero with positive and negative slopes
def setup_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=103, ncols=102, nframes=1, grouptime=1.0,gain=1, deltatime=1):

    times = np.array(list(range(ngroups)),dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    pixdq = np.zeros(shape=(nrows, ncols), dtype= np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
    model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name='MIRI'
    model1.meta.instrument.detector='MIRIMAGE'
    model1.meta.instrument.filter='F480M'
    model1.meta.observation.date='2015-10-13'
    model1.meta.exposure.type='MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name='FULL'
    model1.meta.subarray.xstart=1
    model1.meta.subarray.ystart = 1
    model1.meta.subarray.xsize = ncols
    model1.meta.subarray.ysize = nrows
    model1.meta.exposure.frame_time =deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name='MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = ncols
    gain.meta.subarray.ysize = nrows
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name='MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = ncols
    rnModel.meta.subarray.ysize = nrows
    return model1, gdq, rnModel, pixdq, err, gain


def setup_subarray_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=1032, ncols=1024, subxstart=1, subystart=1,
                 subxsize=1024, subysize=1032, nframes=1,
                 grouptime=1.0,gain=1, deltatime=1):

    times = np.array(list(range(ngroups)),dtype=np.float64) * deltatime
    gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
    err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
    data = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.float64)
    pixdq = np.zeros(shape=(subysize, subxsize), dtype= np.float64)
    read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
    gdq = np.zeros(shape=(nints, ngroups, subysize, subxsize), dtype=np.int32)
    model1 = RampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
    model1.meta.instrument.name='MIRI'
    model1.meta.instrument.detector='MIRIMAGE'
    model1.meta.instrument.filter='F480M'
    model1.meta.observation.date='2015-10-13'
    model1.meta.exposure.type='MIR_IMAGE'
    model1.meta.exposure.group_time = deltatime
    model1.meta.subarray.name='FULL'
    model1.meta.subarray.xstart=subxstart
    model1.meta.subarray.ystart = subystart
    model1.meta.subarray.xsize = subxsize
    model1.meta.subarray.ysize = subysize
    model1.meta.exposure.frame_time =deltatime
    model1.meta.exposure.ngroups = ngroups
    model1.meta.exposure.group_time = deltatime
    model1.meta.exposure.nframes = 1
    model1.meta.exposure.groupgap = 0
    gain = GainModel(data=gain)
    gain.meta.instrument.name='MIRI'
    gain.meta.subarray.xstart = 1
    gain.meta.subarray.ystart = 1
    gain.meta.subarray.xsize = 1024
    gain.meta.subarray.ysize = 1032
    rnModel = ReadnoiseModel(data=read_noise)
    rnModel.meta.instrument.name='MIRI'
    rnModel.meta.subarray.xstart = 1
    rnModel.meta.subarray.ystart = 1
    rnModel.meta.subarray.xsize = 1024
    rnModel.meta.subarray.ysize = 1032
    return model1, gdq, rnModel, pixdq, err, gain
