import gwcs
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from jwst.assign_wcs.niriss import imaging
from jwst.assign_wcs.util import (
    apply_mosaic_member_tweak,
    get_mosaic_member_tweak,
    wfss_imaging_wcs,
)
from jwst.tests.niriss_rate_helpers import niriss_image_rate_model, niriss_wfss_rate_model
from jwst.tests.wcs_helpers import get_reference_files


@pytest.fixture(scope="module")
def tweak_table_rows():
    # Add 3 rows with different visit ids
    rows = []
    for i in range(1, 4):
        rows.append((f"test{i}.fits", str(i), 45.0 + i, 45.0 + i, i, -i, 1, -1 / i, 1 / i, 1))

    # Add one more row with same visit id, different RA/Dec to test multiple matches
    rows.append(("test4.fits", 3, 45.0, 45.0, 4, -4, 1, -1 / 4, 1 / 4, 1))

    # And one more with the identity matrix, for testing shifts
    rows.append(("test5.fits", 5, 45.0, 45.0, 5, -5, 1, 0, 0, 1))
    return rows


@pytest.fixture(scope="module")
def direct_image(tmp_path_factory, tweak_table_rows):
    # Direct image model
    image = niriss_image_rate_model(with_wcs=False)

    # Add tweaks to the image model
    image.tweak_table = np.array(tweak_table_rows, dtype=image.get_dtype("tweak_table"))

    image_path = str(tmp_path_factory.mktemp("tweak_table") / "test_direct.fits")
    image.save(image_path)
    image.close()

    return image_path


@pytest.fixture(scope="module")
def direct_image_no_tweak(tmp_path_factory):
    # Direct image model
    image = niriss_image_rate_model(with_wcs=False)

    image_path = str(tmp_path_factory.mktemp("tweak_table") / "test_direct_no_tweak.fits")
    image.save(image_path)
    image.close()

    return image_path


@pytest.fixture(scope="module")
def wfss_model():
    return niriss_wfss_rate_model(with_wcs=False)


@pytest.fixture(scope="module")
def image_pipeline(wfss_model):
    ref_files = get_reference_files(wfss_model)
    pipe = imaging(wfss_model, ref_files)
    return pipe


@pytest.mark.parametrize("visit_id,expected_row", [("1", 0), ("2", 1), ("3", 3), ("4", None)])
def test_get_mosaic_member_tweak(
    caplog, direct_image, wfss_model, tweak_table_rows, visit_id, expected_row
):
    caplog.set_level("INFO", "jwst")

    # WFSS model to match
    wfss = wfss_model.copy()
    wfss.meta.observation.visit_id = visit_id
    wfss.meta.direct_image = direct_image
    wfss.meta.wcsinfo.ra_ref = 45.0
    wfss.meta.wcsinfo.dec_ref = 45.0

    matrix, shift = get_mosaic_member_tweak(wfss)
    if expected_row is None:
        assert matrix is None
        assert shift is None
        assert "No matching tweak" in caplog.text
    else:
        expected_shift = tweak_table_rows[expected_row][4:6]
        expected_matrix = [
            tweak_table_rows[expected_row][6:8],
            tweak_table_rows[expected_row][8:10],
        ]
        assert np.allclose(matrix, expected_matrix)
        assert np.allclose(shift, expected_shift)
        assert f"Applying WCS tweak from mosaic member test{expected_row + 1}.fits" in caplog.text


@pytest.mark.parametrize("missing_file", [None, "bad.fits"])
def test_get_mosaic_member_tweak_no_file(caplog, wfss_model, missing_file):
    caplog.set_level("INFO", "jwst")

    # No direct image available
    wfss = wfss_model.copy()
    wfss.meta.direct_image = missing_file

    matrix, shift = get_mosaic_member_tweak(wfss)
    assert matrix is None
    assert shift is None
    if missing_file is None:
        assert "No direct image file identified" in caplog.text
    else:
        assert f"Direct image file {missing_file} not found" in caplog.text


def test_get_mosaic_member_tweak_no_table(caplog, wfss_model, direct_image_no_tweak):
    caplog.set_level("INFO", "jwst")

    # No tweak table in direct image
    wfss = wfss_model.copy()
    wfss.meta.direct_image = direct_image_no_tweak

    matrix, shift = get_mosaic_member_tweak(wfss)
    assert matrix is None
    assert shift is None
    assert "Direct image does not contain WCS tweaks" in caplog.text


def test_apply_mosaic_member_tweak(wfss_model, direct_image, tweak_table_rows, image_pipeline):
    # WFSS model to match tweak table
    wfss = wfss_model.copy()
    wfss.meta.observation.visit_id = "1"
    wfss.meta.direct_image = direct_image

    new_pipe = apply_mosaic_member_tweak(wfss, image_pipeline)
    new_wcs = gwcs.WCS(new_pipe)
    assert "v2v3corr" in new_wcs.available_frames

    # matrix and shift from tweak table are stored in tp_affine part of transform
    expected_shift = tweak_table_rows[0][4:6]
    expected_matrix = [tweak_table_rows[0][6:8], tweak_table_rows[0][8:10]]
    rad2arcsec = 3600.0 * np.rad2deg(1.0)
    tfm = new_wcs.get_transform("v2v3", "v2v3corr")
    assert np.allclose(tfm["tp_affine"].matrix.value, expected_matrix)
    assert np.allclose(tfm["tp_affine"].translation.value * rad2arcsec, expected_shift)


def test_apply_mosaic_member_tweak_none_found(caplog, wfss_model, direct_image, image_pipeline):
    caplog.set_level("INFO", "jwst")

    # WFSS model: visit does not match tweak table in image
    wfss = wfss_model.copy()
    wfss.meta.observation.visit_id = "6"
    wfss.meta.direct_image = direct_image

    new_pipe = apply_mosaic_member_tweak(wfss, image_pipeline)
    assert new_pipe is image_pipeline
    new_wcs = gwcs.WCS(new_pipe)
    assert "v2v3corr" not in new_wcs.available_frames

    assert "No matching tweak" in caplog.text


def test_wfss_imaging_wcs_with_tweak(wfss_model, direct_image, tweak_table_rows):
    # wfss model that will pick up a tweak from the direct image
    wfss_tweak = wfss_model.copy()
    wfss_tweak.meta.observation.visit_id = "5"
    wfss_tweak.meta.direct_image = direct_image

    # wfss model that will get no tweak
    wfss_no_tweak = wfss_model.copy()

    # Update fits wcsinfo in both
    wfss_imaging_wcs(wfss_tweak, imaging, max_pix_error=0.05, bbox=((1, 1024), (1, 1024)))
    wfss_imaging_wcs(wfss_no_tweak, imaging, max_pix_error=0.05, bbox=((1, 1024), (1, 1024)))

    for key in ["a_order", "b_order", "crpix1", "crpix2", "crval1", "crval2", "cd1_1"]:
        assert key in wfss_tweak.meta.wcsinfo.instance
        assert key in wfss_no_tweak.meta.wcsinfo.instance

    # Updated position in the tweaked file is at the expected separation from the untweaked file
    coord1 = SkyCoord(wfss_tweak.meta.wcsinfo.crval1, wfss_tweak.meta.wcsinfo.crval2, unit="deg")
    coord2 = SkyCoord(
        wfss_no_tweak.meta.wcsinfo.crval1, wfss_no_tweak.meta.wcsinfo.crval2, unit="deg"
    )
    separation = coord1.separation(coord2).to("arcsec").value
    expected_shift = np.array(tweak_table_rows[-1][4:6])
    expected_sep = np.sqrt(np.sum(expected_shift**2))
    assert np.isclose(separation, expected_sep)
