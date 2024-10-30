import pytest
import numpy as np
from jwst.extract_1d.soss_extract import atoca


@pytest.fixture(scope="module")
def wave_map():
    shp = (100, 30)
    wave_ord1 = np.linspace(2.8, 0.8, shp[0])

    pass

@pytest.fixture(scope="module")
def trace_profile(wave_map):
    pass

@pytest.fixture(scope="module")
def throughput(wave_map):
    pass

@pytest.fixture(scope="module")
def kernels():
    pass

@pytest.fixture(scope="module")
def wave_grid():
    pass

@pytest.fixture(scope="module")
def mask_trace_profile(wave_map):
    pass


def test_extraction_engine():
    pass
