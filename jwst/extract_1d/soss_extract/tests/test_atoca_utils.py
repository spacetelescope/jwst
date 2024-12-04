
import pytest
from jwst.extract_1d.soss_extract import atoca_utils as au
import numpy as np




# wavelengths have min, max (1.5, 4.0) and a bit of non-linearity
WAVELENGTHS = np.linspace(1.5, 3.0, 50) + np.sin(np.linspace(0, np.pi/2, 50))
@pytest.fixture(scope="function")
def wave_map():
    wave_map = np.array([
        WAVELENGTHS,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.2,
        WAVELENGTHS+0.4,
    ])
    wave_map[1,3] = -1 #test skip of bad value
    return wave_map


@pytest.fixture(scope="function")
def wave_map_o2(wave_map):
    return np.copy(wave_map) - 1.0
    

@pytest.fixture(scope="function")
def trace_profile(wave_map):
    thrpt = np.array([0.01, 0.95, 1.0, 0.8, 0.01])
    trace_profile = np.ones_like(wave_map)
    return trace_profile*thrpt[:,None]


@pytest.fixture(scope="function")
def trace_profile_o2(wave_map_o2):
    thrpt = np.array([0.001, 0.01, 0.01, 0.2, 0.99])
    trace_profile = np.ones_like(wave_map_o2)
    return trace_profile*thrpt[:,None]









@pytest.fixture(scope="module")
def kern_array(kernel):

    return au._fct_to_array(kernel, np.linspace(2.0, 4.0, 50), [0, 50], 1e-5)


def test_sparse_c(kern_array):
    """Here kernel must be a 2-D array already, of shape (N_ker, N_k_convolved)"""

    # test typical case n_k = n_kc and i=0
    n_k = kern_array.shape[1]
    i_zero = 0
    matrix = au._sparse_c(kern_array, n_k, i_zero)

    # TODO: add more here


def test_get_c_matrix(kernel):
    """See also test_fct_to_array and test_sparse_c for more detailed tests
    of functions called by this one"""


    #TODO: what to use for grid? is it / can it be the same as wave_trace?
    # I think it can be the same but does not need to be, and would be a better
    # test if it were different, because the kernel is an interpolator that was
    # created using wavelengths that are included in wave_trace.
    matrix = au.get_c_matrix(kernel, grid, i_bounds=None, thresh=1e-5)

    # test with WebbKernel as the kernel
    # ensure normalized
    # ensure sparse


    # test where input kernel is a 2-D array instead of callable


    # test where input kernel is size 1


    # test where i_bounds is not None


    # Test invalid kernel input (wrong dimensions)



@pytest.fixture(scope="module")
def tikhoTests():
    """Make a TikhoTests dictionary"""


    return au.TikhoTests({'factors': factors,
                          'solution': sln,
                          'error': err,
                          'reg': reg,
                          'grid': wave_grid})



def test_tikho_tests(tikhoTests):

    assert False