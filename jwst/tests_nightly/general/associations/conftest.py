"""Pytest configurations"""
import pytest

# Add option to specify a single pool name
def pytest_addoption(parser):
    parser.addoption(
        '--sdp-pool', metavar='sdp_pool', default=None,
        help='SDP pool to run. Specify the name only, not extension or path'
    )


@pytest.fixture
def sdp_pool(request):
    return request.config.getoption('--sdp-pool')
