import pytest
from astropy.modeling.mappings import Identity
from stdatamodels.jwst.datamodels import ChromCorrModel

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep


@pytest.fixture(autouse=True)
def mock_chromcorr_ref(monkeypatch, tmp_cwd):
    """
    Return an identity mock in place of the chromcorr reference file.

    Patches several helper functions that call ``step.get_reference_file`` directly.
    Can be removed in Build 13.2 once the chromcorr file is always found in CRDS.
    """
    _orig = AssignWcsStep.get_reference_file

    def mock_chromcorr(pathname):
        """Make a mock chromcorr file with an identity transform."""  # numpydoc ignore: RT01
        chrom_corr_model = Identity(3)
        chromcorr = ChromCorrModel()
        chromcorr.model = chrom_corr_model
        chromcorr.meta.description = "Mock chromcorr file for testing"
        chromcorr.meta.author = "test"
        chromcorr.meta.pedigree = "GROUND"
        chromcorr.meta.useafter = "2000-01-01T00:00:00"
        fname = pathname / "mock_chromcorr.asdf"
        chromcorr.save(fname)
        return str(fname)

    def _patched(self, model, reference_file_type):
        if reference_file_type == "chromcorr":
            return mock_chromcorr(tmp_cwd)
        return _orig(self, model, reference_file_type)

    monkeypatch.setattr(AssignWcsStep, "get_reference_file", _patched)
