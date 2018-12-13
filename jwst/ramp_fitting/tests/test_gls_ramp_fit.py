import numpy as np

from jwst.ramp_fitting.tests.check_gls_stats import simple_ramp


#ramp_fit_step hardcodes the input to be OLS. So you can't get to the GLS code.
#@pytest.mark.xfail(reason="Fails, not implemented")
def test_gls_ramp_40():
    #Here given a 10 group ramp with an exact slope of 20/group.
    # The output slope should be 20.
    slopes = simple_ramp(fit_method='GLS', ngroups=40, cr_group=None)

    # print(slopes[0].data)
    # take the ratio of the slopes to get the relative error
    # assert_almost_equals(slopes[0].data[500, 500]/ 20.0, 1.0, places=6)
    np.testing.assert_allclose(slopes[0].data[500, 500]/ 20.0, 1.0)
