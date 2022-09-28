"""
Tests of the migrate_data script, which attempts to update
files that have become invalid due to changes in model schemas
between versions of this package.

Obtains examples of files from artifactory truth data for
older releases.
"""
import subprocess

from astropy.io import fits
import pytest

from jwst import datamodels


@pytest.fixture(autouse=True)
def strict_validation(monkeypatch):
    monkeypatch.setenv("STRICT_VALIDATION", "true")
    yield


@pytest.mark.skip(reason="obsolete test with no available old input data")
@pytest.mark.bigdata
@pytest.mark.parametrize("truth_path", [
    "truth/test_miri_lrs_slit_spec2/jw00623032001_03102_00001_mirimage_x1d.fits",
    "truth/test_nirspec_mos_spec2/f170lp-g235m_mos_observation-6-c0e0_001_dn_nrs1_mod_x1d.fits",
], ids=["miri-lrs-x1d", "nirspec-mos-x1d"])
def test_x1d_spec_table(truth_path, rtdata):
    rtdata.env = "1.1.0"
    rtdata.get_truth(truth_path)

    # Confirm that the file doesn't initially validate
    # (open with fits first so that the failed call to open doesn't leave behind an open file)
    with fits.open(rtdata.truth, memmap=False) as hdu:
        with pytest.raises(ValueError, match="Column names don't match schema"):
            with datamodels.open(hdu):
                pass

    subprocess.check_call(["migrate_data", rtdata.truth, "--in-place"])

    # Now the model should validate
    with datamodels.open(rtdata.truth):
        pass


@pytest.mark.skip(reason="obsolete test with no available old input data")
@pytest.mark.bigdata
@pytest.mark.parametrize("truth_path",
                         ["truth/test_nircam_mtimage/mt_pre_1-2-2_schema_uncal.fits"],
                         ids=["nircam_mt"])
def test_pre_1_2_2_schemas(truth_path, rtdata):
    rtdata.env = "1.2.1"
    rtdata.get_truth(truth_path)

    # Confirm that the file doesn't initially validate
    # (open with fits first so that the failed call to open doesn't leave behind an open file)
    with fits.open(rtdata.truth, memmap=False) as hdu:
        with pytest.raises(ValueError, match="Column names don't match schema"):
            with datamodels.open(hdu):
                pass

    subprocess.check_call(["migrate_data", rtdata.truth, "--in-place"])

    # Now the model should validate
    with datamodels.open(rtdata.truth):
        pass
