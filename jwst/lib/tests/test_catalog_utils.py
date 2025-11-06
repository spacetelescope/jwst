import pytest
from astropy.utils.diff import report_diff_values

from jwst.lib.catalog_utils import read_source_catalog
from jwst.source_catalog import SourceCatalogStep
from jwst.source_catalog.tests.helpers import make_nircam_model


@pytest.fixture
def nircam_model():
    return make_nircam_model()


def test_read_source_catalog(nircam_model, tmp_path):
    """Ensure read_source_catalog can read a real source catalog file."""

    step = SourceCatalogStep(
        snr_threshold=0.5, npixels=10, bkg_boxsize=50, kernel_fwhm=2.0, save_results=False
    )
    cat = step.run(nircam_model)

    # in memory
    cat_in_memory = read_source_catalog(cat)
    assert len(cat_in_memory) > 0

    # from file
    catpath = tmp_path / "test.ecsv"
    catstr = str(catpath)
    cat.write(catpath, format="ascii.ecsv", overwrite=True)
    cat_from_str = read_source_catalog(catstr)
    identical = report_diff_values(cat_in_memory, cat_from_str)
    assert identical, f"Catalogs read in memory and from file str differ: {identical}"


def test_read_source_catalog_invalid_input():
    """Test that read_source_catalog raises appropriate errors for invalid inputs."""

    # Test empty string
    with pytest.raises(ValueError, match="Empty catalog filename"):
        read_source_catalog("")

    # Test file not found
    with pytest.raises(FileNotFoundError, match="Could not find catalog"):
        read_source_catalog("nonexistent_catalog.ecsv")

    # Test invalid type
    with pytest.raises(
        TypeError, match="Need to input string name of catalog or astropy.table.table.QTable"
    ):
        read_source_catalog(12345)
