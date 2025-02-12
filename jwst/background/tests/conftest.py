import pytest
import pathlib


@pytest.fixture(scope="module")
def data_path():
    return pathlib.Path(__file__).parent / "data"
