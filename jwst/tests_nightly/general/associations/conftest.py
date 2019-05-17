"""Pytest configurations"""
import pytest

from jwst.tests_nightly.general.associations.sdp_pools_source import SDPPoolsSource


# Add option to specify a single pool name
def pytest_addoption(parser):
    parser.addoption(
        '--sdp-pool', metavar='sdp_pool', default=None,
        help='SDP pool to run. Specify the name only, not extension or path'
    )


@pytest.fixture
def sdp_pool(request):
    return request.config.getoption('--sdp-pool')


def pytest_generate_tests(metafunc):
    """Prefetch and parametrize a set of test pools"""
    if 'pool_path' in metafunc.fixturenames:
        SDPPoolsSource.inputs_root = metafunc.config.getini('inputs_root')[0]
        SDPPoolsSource.results_root = metafunc.config.getini('results_root')[0]

        pools = SDPPoolsSource()

        try:
            pool_paths = pools.pool_paths
        except Exception:
            pool_paths = []

        metafunc.parametrize('pool_path', pool_paths)
