import pytest
import numpy as np

from jwst.ramp_fitting.tests.check_gls_stats import simple_ramp


@pytest.mark.parametrize("method", ['GLS', 'OLS'])
@pytest.mark.parametrize("ngroups", [10, 20])
@pytest.mark.parametrize("cr_group", [None])
def test_gls_ramp_10(method, ngroups, cr_group):
    slopes = simple_ramp(fit_method=method, ngroups=ngroups,
                         cr_group=cr_group)

    np.testing.assert_allclose(slopes[0].data[500, 500]/ 20.0, 1.0)
