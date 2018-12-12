# program to check GLS cpu and memory usage

import time

from jwst.datamodels import dqflags

from jwst.ramp_fitting.ramp_fit import ramp_fit

from jwst.ramp_fitting.tests.test_gls_ramp_fit import setup_inputs


def simple_ramp(fit_method='OLS'):
    #Here given a 10 group ramp with an exact slope of 20/group.
    # The output slope should be 20.
    ngroups = 5
    slope_per_group = 20.
    model1, gdq, rnModel, pixdq, err, gain = setup_inputs(ngroups=ngroups)
    for k in range(ngroups):
        model1.data[0, k, :, :] = 10.0 + k*slope_per_group

    model1.groupdq[0, 2, :, :] = dqflags.group['JUMP_DET']
    slopes = ramp_fit(model1, 1024*10., True, rnModel, gain, fit_method,
                      'optimal')
    return slopes


if __name__ == '__main__':
    start_time = time.process_time()
    ols_slopes = simple_ramp(fit_method='OLS')
    end_time_ols = time.process_time()
    gls_slopes = simple_ramp(fit_method='GLS')
    end_time_gls = time.process_time()
    print('OLS:', end_time_ols - start_time)
    print('GLS:', end_time_gls - end_time_ols)
