import numpy as np
import pytest
from stdatamodels.jwst import datamodels

from jwst.lib.reffile_utils import MatchRowError
from jwst.pfpc import pfpc as pf


def test_find_correction_no_match(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    # input has channel = '12' prior to cube_build,
    # so there is no matching row in the table
    row = pf.find_correction(mrs_dith1_ch12_medium, mrs_pfpc_model.pfpc_table)
    assert row is None
    assert "Expected to find one matching row in table, found 0" in caplog.text


@pytest.mark.parametrize("require_one", [True, False])
def test_find_correction_require_one(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model, require_one):
    # if we ignore channel, there will be multiple matches in the table
    if require_one:
        # raises an error
        with pytest.raises(
            MatchRowError, match="Expected to find one matching row in table, found 4"
        ):
            pf.find_correction(
                mrs_dith1_ch12_medium,
                mrs_pfpc_model.pfpc_table,
                ignore_keys=["channel"],
                require_one=require_one,
            )

    else:
        # returns the first matching row
        row = pf.find_correction(
            mrs_dith1_ch12_medium,
            mrs_pfpc_model.pfpc_table,
            ignore_keys=["channel"],
            require_one=require_one,
        )
        assert len(row["channel"]) == 1


def test_find_correction_match(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    # set input channel to 2, so it finds one matching row
    input_copy = mrs_dith1_ch12_medium.copy()
    input_copy.meta.instrument.channel = "2"

    row = pf.find_correction(input_copy, mrs_pfpc_model.pfpc_table)
    assert row["channel"] == "2"


def test_find_correction_match_reference(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    input_copy = mrs_dith1_ch12_medium.copy()
    input_copy.meta.instrument.channel = "2"
    input_copy.meta.ref_file.photom.name = "crds://test_photom_1.fits"

    # Table value matches input
    pfpc_copy = mrs_pfpc_model.copy()
    pfpc_copy.pfpc_table["r_photom"] = "test_photom_1.fits"

    # A matching row is returned
    row = pf.find_correction(input_copy, mrs_pfpc_model.pfpc_table)
    assert row["channel"] == "2"

    # No warning: the file name matches
    assert "does not match input value" not in caplog.text


def test_find_correction_mismatched_reference(caplog, mrs_dith1_ch12_medium, mrs_pfpc_model):
    input_copy = mrs_dith1_ch12_medium.copy()
    input_copy.meta.instrument.channel = "2"
    input_copy.meta.ref_file.photom.name = "crds://test_photom_1.fits"

    # Table value does not match input
    pfpc_copy = mrs_pfpc_model.copy()
    pfpc_copy.pfpc_table["r_photom"] = "test_photom_2.fits"

    # A matching row is still returned
    row = pf.find_correction(input_copy, mrs_pfpc_model.pfpc_table)
    assert row["channel"] == "2"

    # Warns: the file name does not match
    assert "test_photom_2.fits does not match input value test_photom_1.fits" in caplog.text
    assert "PFPC corrections may be invalid" in caplog.text


def test_process_unsupported_exptype():
    model = datamodels.ImageModel()
    model.meta.exposure.type = "NRC_IMAGE"
    with pytest.raises(TypeError, match="Exposure type NRC_IMAGE is not supported"):
        pf.process_exposures([model])


def test_apply_correction(caplog, mrs_dith1_ch1_medium_x1d, mrs_pfpc_model):
    spec_by_exposure = [mrs_dith1_ch1_medium_x1d.copy()]
    pfpc_table = mrs_pfpc_model.pfpc_table
    corrected = pf.apply_correction(spec_by_exposure, pfpc_table)

    assert len(corrected) == 1
    assert len(corrected[0].spec) == 1
    input_spec = mrs_dith1_ch1_medium_x1d.spec[0].spec_table
    corrected_spec = corrected[0].spec[0].spec_table

    # wavelength should be unmodified
    assert np.allclose(corrected_spec["WAVELENGTH"], input_spec["WAVELENGTH"])

    # expected correction in mock table is channel + dither = 2.0
    correction = 2.0
    assert np.allclose(corrected_spec["FLUX"], input_spec["FLUX"] / correction)
    assert np.allclose(corrected_spec["FLUX_ERROR"], input_spec["FLUX_ERROR"] / correction)
    assert np.allclose(corrected_spec["SURF_BRIGHT"], input_spec["SURF_BRIGHT"] / correction)
    assert np.allclose(corrected_spec["SB_ERROR"], input_spec["SB_ERROR"] / correction)
