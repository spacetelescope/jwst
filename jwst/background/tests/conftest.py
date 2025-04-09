import pytest
import pathlib


@pytest.fixture(scope="module")
def data_path():
    """
    Define the data path for testing.

    Returns
    -------
    pathlib object
        The path to where the data lives
    """
    return pathlib.Path(__file__).parent / "data"
