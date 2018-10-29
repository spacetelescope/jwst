import pytest
import numpy as np

from jwst.ramp_fitting.ramp_fit import ramp_fit
from jwst.datamodels import dqflags
from jwst.datamodels import MIRIRampModel
from jwst.datamodels import GainModel, ReadnoiseModel
from nose.tools import assert_almost_equals

#@pytest.mark.skip(reason="not using now")
def test_nocrs_noflux():
        #all pixel values are zero. So slope should be zero
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5)
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    print(slopes[0].data)
    assert(0 == np.max(slopes[0].data))
    assert(0 == np.min(slopes[0].data))
    
#@pytest.mark.skip(reason="not using now")
def test_simple_ramp():
    #Here given a 10 group ramp with an exact slope of 20/group. The output slope should be 20.
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=10, deltatime=3)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 30.0
    model1.data[0, 2, 500, 500] = 50.0
    model1.data[0, 3, 500, 500] = 70.0
    model1.data[0, 4, 500, 500] = 90.0
    model1.data[0, 5, 500, 500] = 110.0
    model1.data[0, 6, 500, 500] = 130.0
    model1.data[0, 7, 500, 500] = 150.0
    model1.data[0, 8, 500, 500] = 170.0
    model1.data[0, 9, 500, 500] = 190.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    print(slopes[0].data[500,500])
    # take the ratio of the slopes to get the relative error
    assert_almost_equals(slopes[0].data[500, 500]/ (20.0/3), 1.0, places=6)

#@pytest.mark.skip(reason="not using now")
def test_read_noise_only_fit():
        #
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,readnoise=50)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 15.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 33.0
    model1.data[0, 4, 500, 500] = 60.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    xvalues = np.arange(5)*1.0
    yvalues = np.array([10,15,25,33,60])
    coeff = np.polyfit(xvalues, yvalues, 1)
    print('x ',repr(xvalues))
    print('coeff ',repr(coeff))
    assert_almost_equals(slopes[0].data[500, 500]/coeff[0], 1.0, places=6)

#@pytest.mark.skip(reason="not using now")
def test_photon_noise_only_fit():
        #
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=5,gain=1000,readnoise=1)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 15.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 33.0
    model1.data[0, 4, 500, 500] = 60.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    cds_slope = (model1.data[0,4,500,500] - model1.data[0,0,500,500])/ 4.0
    print('CDS Slope', cds_slope)
    print('Slope ',slopes[0].data[500, 500])

    assert_almost_equals(slopes[0].data[500, 500]/cds_slope, 1.0, places=2)

#@pytest.mark.skip(reason="not using now")
def test_two_groups_fit():
        #
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=2,gain=1,readnoise=10)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 15.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    cds_slope = (model1.data[0,1,500,500] - model1.data[0,0,500,500])
    print('CDS Slope', cds_slope)
    print('Slope ',slopes[0].data[500, 500])

    assert_almost_equals(slopes[0].data[500, 500]/cds_slope, 1.0, places=2)

#ramp_fit_step hardcodes the input to be OLS. So you can't get to the GLS code.
@pytest.mark.xfail(reason="Fails, not implemented")
def test_simple_gls_ramp():
    #Here given a 10 group ramp with an exact slope of 20/group. The output slope should be 20.
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=10)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 30.0
    model1.data[0, 2, 500, 500] = 50.0
    model1.data[0, 3, 500, 500] = 70.0
    model1.data[0, 4, 500, 500] = 90.0
    model1.data[0, 5, 500, 500] = 110.0
    model1.data[0, 6, 500, 500] = 130.0
    model1.data[0, 7, 500, 500] = 150.0
    model1.data[0, 8, 500, 500] = 170.0
    model1.data[0, 9, 500, 500] = 190.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'GLS', 'optimal')
    print(slopes[0].data)
    # take the ratio of the slopes to get the relative error
    assert_almost_equals(slopes[0].data[500, 500]/ 20.0, 1.0, places=6)

#@pytest.mark.skip(reason="not using now")
def test_two_groups_unc():
    grouptime=3.0
    deltaDN = 5
    ingain = 2
    inreadnoise =10
    ngroups=2
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups, gain=ingain, readnoise=inreadnoise,deltatime=grouptime)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 10.0 + deltaDN
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    print('Slope ',slopes[0].data[500, 500])
    print(' Poisson Variance ',slopes[0].var_poisson[500,500])
    print(' Read Variance ',slopes[0].var_rnoise[500,500])
    print(' Total Sigma ',slopes[0].err[500,500])
    delta_electrons = deltaDN * ingain
    single_sample_readnoise = inreadnoise/np.sqrt(2)
    print('CDS/ SS Readnoise',inreadnoise/single_sample_readnoise)
    print('CDS variance',(inreadnoise**2/grouptime**2))
    print('Slope Variance',(12 * single_sample_readnoise**2)/(ngroups*(ngroups**2 - 1)*grouptime**2))
    
    assert_almost_equals(slopes[0].var_poisson[500,500]/((deltaDN/ingain)/grouptime**2), 1.0, places=7)
    assert_almost_equals(slopes[0].var_rnoise[500,500]/(inreadnoise**2/grouptime**2), 1.0, places=7)
    assert_almost_equals(slopes[0].var_rnoise[500,500]/(12*single_sample_readnoise**2/(ngroups*(ngroups**2 - 1)*grouptime**2)), 1.0, places=7)
    assert_almost_equals(slopes[0].err[500,500]/(np.sqrt((deltaDN/ingain)/grouptime**2+(inreadnoise**2/grouptime**2))), 1.0, places=7)

#@pytest.mark.skip(reason="not using now")
def test_five_groups_unc():
    grouptime=3.0
    deltaDN = 5
    ingain = 2
    inreadnoise =7
    ngroups=5
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise,deltatime=grouptime)
    model1.data[0, 0, 500, 500] = 10.0
    model1.data[0, 1, 500, 500] = 15.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 33.0
    model1.data[0, 4, 500, 500] = 60.0
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, 'OLS', 'optimal')
    print('Slope ',slopes[0].data[500, 500])
    print(' Poisson Variance ',slopes[0].var_poisson[500,500])
    print(' Read Variance ',slopes[0].var_rnoise[500,500])
    print(' Total Sigma ',slopes[0].err[500,500])
    out_slope=slopes[0].data[500, 500]
    median_slope=np.median(np.diff(model1.data[0,:,500,500]))/grouptime
    print('median slope',median_slope)
    deltaDN = 50
    delta_time = (ngroups - 1) * grouptime
    delta_electrons = median_slope * ingain *delta_time
    print('delta electons',delta_electrons)
    print('delta time', delta_time)
    single_sample_readnoise = np.float64(inreadnoise/np.sqrt(2))
    assert_almost_equals(slopes[0].var_poisson[500,500]/((median_slope)/(ingain*delta_time)), 1.0, places=6)
    assert_almost_equals(slopes[0].var_rnoise[500,500]/(12 * single_sample_readnoise**2/(ngroups * (ngroups**2 - 1) * grouptime**2)), 1.0, places=7)
    assert_almost_equals(slopes[0].err[500,500]/np.sqrt(slopes[0].var_poisson[500,500]  + slopes[0].var_rnoise[500,500] ), 1.0, places=6)

##@pytest.mark.skip(reason="not using now")
def test_oneCR_10_groups_combination():
    grouptime=3.0
    deltaDN = 5
    ingain = 200 # use large gain to show that Poisson noise doesn't affect the recombination
    inreadnoise = np.float64(7)
    ngroups=10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise,deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 500, 500] = 15.0
    model1.data[0, 1, 500, 500] = 20.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 30.0
    model1.data[0, 4, 500, 500] = 35.0
    model1.data[0, 5, 500, 500] = 140.0
    model1.data[0, 6, 500, 500] = 150.0
    model1.data[0, 7, 500, 500] = 160.0
    model1.data[0, 8, 500, 500] = 170.0
    model1.data[0, 9, 500, 500] = 180.0
    model1.groupdq[0,5,500,500]=dqflags.group['JUMP_DET']
    slopes, int_model, opt_model, gls_opt_model= ramp_fit(model1, 64000,  True, rnModel, gain, 'OLS', 'optimal')
    print("slopes shape",type(slopes))
    print("opt_model shape",type(opt_model))
    print("int_model shape",type(int_model))
    print("gls_model shape",type(gls_opt_model))
    print('Slope ',slopes.data[500, 500])
    print(' Poisson Variance ',slopes.var_poisson[500,500])
    print(' Read Variance ',slopes.var_rnoise[500,500])
    print(' Total Sigma ',slopes.err[500,500])
    print(' int model   exposure slope',int_model.data[0,500,500])
    print(' int model  seg 0 slope',opt_model.slope[0,0,500,500])
    print(' int model  seg 1 slope',opt_model.slope[0,1,500,500])
    print(' int model   seg 0 RN var',opt_model.var_rnoise[0,0,500,500])
    print(' int model   seg 1 RN var',opt_model.var_rnoise[0,1,500,500])
    print(' int model   crmag',opt_model.crmag[0,0,500,500])
    print(' int model  time type',type(int_model.int_times))

    segment_groups  = 5
    single_sample_readnoise = np.float64( inreadnoise/np.sqrt(2))
    #check that the segment variance is as expected
    np.testing.assert_allclose(opt_model.var_rnoise[0,0,500,500],(12.0 * single_sample_readnoise**2/(segment_groups * (segment_groups**2 - 1) * grouptime**2)), rtol=1e-6)
    # check the combined slope is the average of the two segments since they have the same number of groups 
    np.testing.assert_allclose(slopes.data[500, 500], 2.5,rtol=1e-5)
    #check that the slopes of the two segments are correct
    np.testing.assert_allclose(opt_model.slope[0,0,500, 500], 5/3.0,rtol=1e-5)
    np.testing.assert_allclose(opt_model.slope[0,1,500, 500], 10/3.0,rtol=1e-5)

#@pytest.mark.skip(reason="not using now")    
def test_oneCR_10_groups_combination_noisy2ndSegment():
    grouptime=3.0
    deltaDN = 5
    ingain = 200 # use large gain to show that Poisson noise doesn't affect the recombination
    inreadnoise =7
    ngroups=10
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups,
                                                          gain=ingain, readnoise=inreadnoise,deltatime=grouptime)
    # two segments perfect fit, second segment has twice the slope
    model1.data[0, 0, 500, 500] = 15.0
    model1.data[0, 1, 500, 500] = 20.0
    model1.data[0, 2, 500, 500] = 25.0
    model1.data[0, 3, 500, 500] = 30.0
    model1.data[0, 4, 500, 500] = 35.0
    model1.data[0, 5, 500, 500] = 135.0
    model1.data[0, 6, 500, 500] = 155.0
    model1.data[0, 7, 500, 500] = 160.0
    model1.data[0, 8, 500, 500] = 168.0
    model1.data[0, 9, 500, 500] = 180.0
    model1.groupdq[0,5,500,500]=dqflags.group['JUMP_DET']
    slopes, int_model, opt_model, gls_opt_model= ramp_fit(model1, 64000,  True, rnModel, gain, 'OLS', 'optimal')
    print("slopes shape",type(slopes))
    print("opt_model shape",type(opt_model))
    print("int_model shape",type(int_model))
    print("gls_model shape",type(gls_opt_model))
    print('Slope ',slopes.data[500, 500])
    print(' Poisson Variance ',slopes.var_poisson[500,500])
    print(' Read Variance ',slopes.var_rnoise[500,500])
    print(' Total Sigma ',slopes.err[500,500])
    print(' int model   exposure slope',int_model.data[0,500,500])
    print(' int model  seg 0 slope',opt_model.slope[0,0,500,500])
    print(' int model  seg 1 slope',opt_model.slope[0,1,500,500])
    print(' int model   seg 0 RN var',opt_model.var_rnoise[0,0,500,500])
    print(' int model   seg 1 RN var',opt_model.var_rnoise[0,1,500,500])
    print(' int model   crmag',opt_model.crmag[0,0,500,500])
    print(' int model  time type',type(int_model.int_times))

    avg_slope = (opt_model.slope[0,0,500, 500] + opt_model.slope[0,1,500, 500])/2.0
        #even with noiser second segment, final slope should be just the average since they have the same number of groups
    np.testing.assert_allclose(slopes.data[500, 500], avg_slope,rtol=1e-5)


    
#Need test for multi-ints near zero with positive and negative slopes

def setup_inputs(ngroups=10, readnoise=10, nints=1,
                 nrows=1032, ncols=1024, nframes=1, grouptime=1.0,gain=1, deltatime=1):
        
        print('readnoise', readnoise)
        print('gain',gain)
        times = np.array(list(range(ngroups)),dtype=np.float64) * deltatime
        gain = np.ones(shape=(nrows, ncols), dtype=np.float64) * gain
        err = np.ones(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        data = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.float64)
        pixdq = np.zeros(shape=(nrows, ncols), dtype= np.float64)
        read_noise = np.full((nrows, ncols), readnoise, dtype=np.float64)
        gdq = np.zeros(shape=(nints, ngroups, nrows, ncols), dtype=np.int32)
        model1 = MIRIRampModel(data=data, err=err, pixeldq=pixdq, groupdq=gdq, times=times)
        model1.meta.instrument.name='MIRI'
        model1.meta.instrument.detector='MIRIMAGE'
        model1.meta.instrument.filter='F480M'
        model1.meta.observation.date='2015-10-13'
        model1.meta.exposure.type='MIR_IMAGE'
        model1.meta.exposure.group_time = deltatime
        model1.meta.subarray.name='FULL'
        model1.meta.subarray.xstart=1
        model1.meta.subarray.ystart = 1
        model1.meta.subarray.xsize = 1024
        model1.meta.subarray.ysize = 1032
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

