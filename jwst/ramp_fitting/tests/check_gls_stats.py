# program to check GLS cpu and memory usage

import time

from jwst.ramp_fitting.ramp_fit import ramp_fit

from jwst.ramp_fitting.tests.test_gls_ramp_fit import setup_inputs


def simple_ramp(fit_method='OLS'):
    #Here given a 10 group ramp with an exact slope of 20/group.
    # The output slope should be 20.
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
    slopes = ramp_fit(model1, 64000, True, rnModel, gain, fit_method,
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
