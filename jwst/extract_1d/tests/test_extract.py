import json

import numpy as np
import pytest
import stdatamodels.jwst.datamodels as dm
from astropy.modeling import polynomial
from numpy.testing import assert_allclose, assert_equal

from jwst.datamodels import ModelContainer
from jwst.extract_1d import extract as ex
from jwst.extract_1d import psf_profile as pp


@pytest.fixture()
def extract1d_ref_dict():
    apertures = [
        {"id": "slit1"},
        {"id": "slit2", "region_type": "other"},
        {"id": "slit3", "xstart": 10, "xstop": 20, "ystart": 10, "ystop": 20},
        {"id": "slit4", "bkg_coeff": [[10], [20]]},
        {"id": "slit5", "bkg_coeff": None},
        {"id": "slit6", "use_source_posn": True},
        {"id": "slit7", "spectral_order": 20},
        {"id": "S200A1"},
        {"id": "S1600A1", "use_source_posn": False},
    ]
    ref_dict = {"apertures": apertures}
    return ref_dict


@pytest.fixture()
def extract1d_ref_file(tmp_path, extract1d_ref_dict):
    filename = str(tmp_path / "extract1d_ref.json")
    with open(filename, "w") as fh:
        json.dump(extract1d_ref_dict, fh)
    return filename


@pytest.fixture()
def extract_defaults():
    default = {
        "bkg_coeff": None,
        "bkg_fit": None,
        "bkg_order": 0,
        "extract_width": None,
        "extraction_type": "box",
        "independent_var": "pixel",
        "match": "exact match",
        "smoothing_length": 0,
        "spectral_order": 1,
        "src_coeff": None,
        "subtract_background": False,
        "position_offset": 0.0,
        "trace": None,
        "use_source_posn": False,
        "model_nod_pair": False,
        "optimize_psf_location": False,
        "xstart": 0,
        "xstop": 49,
        "ystart": 0,
        "ystop": 49,
        "psf": "N/A",
    }
    return default


@pytest.fixture()
def create_extraction_inputs(mock_nirspec_fs_one_slit, extract1d_ref_dict):
    input_model = mock_nirspec_fs_one_slit
    slit = None
    output_model = dm.MultiSpecModel()
    ref_dict = extract1d_ref_dict
    slitname = "S200A1"
    sp_order = 1
    exp_type = "NRS_FIXEDSLIT"
    yield [input_model, slit, output_model, ref_dict, slitname, sp_order, exp_type]
    output_model.close()


def test_read_extract1d_ref(extract1d_ref_dict, extract1d_ref_file):
    ref_dict = ex.read_extract1d_ref(extract1d_ref_file)
    assert ref_dict == extract1d_ref_dict


def test_read_extract1d_ref_bad_json(tmp_path):
    filename = str(tmp_path / "bad_ref.json")
    with open(filename, "w") as fh:
        fh.write("apertures: [bad,]\n")

    with pytest.raises(RuntimeError, match="Invalid JSON extract1d reference"):
        ex.read_extract1d_ref(filename)


def test_read_extract1d_ref_bad_type(tmp_path):
    filename = str(tmp_path / "bad_ref.fits")
    with open(filename, "w") as fh:
        fh.write("bad file\n")

    with pytest.raises(RuntimeError, match="must be JSON"):
        ex.read_extract1d_ref(filename)


def test_read_extract1d_ref_na():
    ref_dict = ex.read_extract1d_ref("N/A")
    assert ref_dict is None


def test_read_apcorr_ref():
    apcorr_model = ex.read_apcorr_ref(None, "MIR_LRS-FIXEDSLIT")
    assert isinstance(apcorr_model, dm.MirLrsApcorrModel)


def test_get_extract_parameters_default(
    mock_nirspec_fs_one_slit, extract1d_ref_dict, extract_defaults
):
    input_model = mock_nirspec_fs_one_slit

    # match a bare entry
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit1", 1, input_model.meta
    )

    # returned value has defaults except that use_source_posn
    # is switched on for NRS_FIXEDSLIT
    expected = extract_defaults
    expected["use_source_posn"] = True

    assert params == expected


def test_get_extract_parameters_na(mock_nirspec_fs_one_slit, extract_defaults):
    input_model = mock_nirspec_fs_one_slit

    # no reference input: defaults returned
    params = ex.get_extract_parameters(None, input_model, "slit1", 1, input_model.meta)
    assert params == extract_defaults


@pytest.mark.parametrize("bgsub", [None, True])
@pytest.mark.parametrize("bgfit", ["poly", "mean", None])
def test_get_extract_parameters_background(
    mock_nirspec_fs_one_slit, extract1d_ref_dict, bgsub, bgfit
):
    input_model = mock_nirspec_fs_one_slit

    # match a slit with background defined
    params = ex.get_extract_parameters(
        extract1d_ref_dict,
        input_model,
        "slit4",
        1,
        input_model.meta,
        subtract_background=bgsub,
        bkg_fit=bgfit,
    )

    # returned value has background switched on
    assert params["subtract_background"] is True
    assert params["bkg_coeff"] is not None
    assert params["bkg_fit"] == "poly"
    assert params["bkg_order"] == 0


def test_get_extract_parameters_bg_ignored(mock_nirspec_fs_one_slit, extract1d_ref_dict):
    input_model = mock_nirspec_fs_one_slit

    # match a slit with background defined
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit4", 1, input_model.meta, subtract_background=False
    )

    # background parameters are ignored
    assert params["subtract_background"] is False
    assert params["bkg_fit"] is None


@pytest.mark.parametrize("slit", ["slit2", "no_match"])
def test_get_extract_parameters_no_match(mock_nirspec_fs_one_slit, extract1d_ref_dict, slit):
    input_model = mock_nirspec_fs_one_slit

    # no slit with an appropriate region_type matched
    params = ex.get_extract_parameters(extract1d_ref_dict, input_model, slit, 1, input_model.meta)
    assert params == {"match": ex.NO_MATCH}


def test_get_extract_parameters_source_posn_exptype(mock_nirspec_bots, extract1d_ref_dict):
    input_model = mock_nirspec_bots
    input_model.meta.exposure.type = "NRS_LAMP"

    # match a bare entry
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit1", 1, input_model.meta, use_source_posn=None
    )

    # use_source_posn is switched off for unknown exptypes
    assert params["use_source_posn"] is False


def test_get_extract_parameters_source_posn_from_ref(mock_nirspec_bots, extract1d_ref_dict):
    input_model = mock_nirspec_bots

    # match an entry that explicitly sets use_source_posn
    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit6", 1, input_model.meta, use_source_posn=None
    )

    # returned value has use_source_posn switched off by default
    # for NRS_BRIGHTOBJ, but ref file overrides
    assert params["use_source_posn"] is True


@pytest.mark.parametrize("length", [3, 4, 2.8, 3.5])
def test_get_extract_parameters_smoothing(mock_nirspec_fs_one_slit, extract1d_ref_dict, length):
    input_model = mock_nirspec_fs_one_slit

    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit1", 1, input_model.meta, smoothing_length=length
    )

    # returned value has input smoothing length, rounded to an
    # odd integer if necessary
    assert params["smoothing_length"] == 3


@pytest.mark.parametrize("length", [-1, 1, 2, 1.3])
def test_get_extract_parameters_smoothing_bad_value(
    mock_nirspec_fs_one_slit, extract1d_ref_dict, length
):
    input_model = mock_nirspec_fs_one_slit

    params = ex.get_extract_parameters(
        extract1d_ref_dict, input_model, "slit1", 1, input_model.meta, smoothing_length=length
    )

    # returned value has smoothing length 0
    assert params["smoothing_length"] == 0


@pytest.mark.parametrize("use_source", [None, True, False])
def test_get_extract_parameters_extraction_type_none(
    mock_nirspec_fs_one_slit, extract1d_ref_dict, use_source, log_watcher
):
    input_model = mock_nirspec_fs_one_slit

    watcher = log_watcher("jwst.extract_1d.extract", message="Using extraction type")
    params = ex.get_extract_parameters(
        extract1d_ref_dict,
        input_model,
        "slit1",
        1,
        input_model.meta,
        extraction_type=None,
        use_source_posn=use_source,
        psf_ref_name="available",
    )
    watcher.assert_seen()

    # Extraction type is set to optimal if use_source_posn is True
    if use_source is None or use_source is True:
        assert params["use_source_posn"] is True
        assert params["extraction_type"] == "optimal"
    else:
        assert params["use_source_posn"] is False
        assert params["extraction_type"] == "box"


@pytest.mark.parametrize("extraction_type", [None, "box", "optimal"])
def test_get_extract_parameters_no_psf(
    mock_nirspec_fs_one_slit, extract1d_ref_dict, extraction_type, log_watcher
):
    input_model = mock_nirspec_fs_one_slit

    watcher = log_watcher("jwst.extract_1d.extract", message="Setting extraction type to 'box'")
    params = ex.get_extract_parameters(
        extract1d_ref_dict,
        input_model,
        "slit1",
        1,
        input_model.meta,
        extraction_type=extraction_type,
        psf_ref_name="N/A",
    )

    # Warning message issued if extraction type was not already 'box'
    if extraction_type != "box":
        watcher.assert_seen()
    else:
        watcher.assert_not_seen()

    # Extraction type is always box if no psf is available
    assert params["extraction_type"] == "box"


def test_log_params(extract_defaults, log_watcher):
    watcher = log_watcher("jwst.extract_1d.extract", message="Extraction parameters")

    # Defaults don't have dispaxis assigned yet - parameters are not logged
    ex.log_initial_parameters(extract_defaults)
    watcher.assert_not_seen()

    # Add dispaxis: parameters are now logged
    extract_defaults["dispaxis"] = 1
    ex.log_initial_parameters(extract_defaults)
    watcher.assert_seen()


def test_create_poly():
    coeff = [1, 2, 3]
    poly = ex.create_poly(coeff)
    assert isinstance(poly, polynomial.Polynomial1D)
    assert poly.degree == 2
    assert poly(2) == 1 + 2 * 2 + 3 * 2**2


def test_create_poly_empty():
    coeff = []
    assert ex.create_poly(coeff) is None


def test_populate_time_keywords(mock_nirspec_bots, mock_10_multi_int_spec):
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_multi_int_spec)

    # ensure time keywords were added to output table
    for i, spec in enumerate(mock_10_multi_int_spec.spec[0].spec_table):
        assert spec["INT_NUM"] == i + 1
        assert spec["MJD-BEG"] == mock_nirspec_bots.int_times["int_start_MJD_UTC"][i]
        assert spec["TDB-END"] == mock_nirspec_bots.int_times["int_end_BJD_TDB"][i]


def test_populate_time_keywords_no_table(
    mock_nirspec_fs_one_slit, mock_10_multi_int_spec, log_watcher
):
    watcher = log_watcher("jwst.extract_1d.extract", message="no INT_TIMES table")
    ex.populate_time_keywords(mock_nirspec_fs_one_slit, mock_10_multi_int_spec)

    # No int_times table: warns and integration is set to simple index
    watcher.assert_seen()
    for spec in mock_10_multi_int_spec.spec:
        assert_equal(spec.spec_table["INT_NUM"], np.arange(1, 11))


def test_populate_time_keywords_multislit(mock_nirspec_mos, mock_10_multi_int_spec):
    mock_nirspec_mos.meta.exposure.nints = 10
    ex.populate_time_keywords(mock_nirspec_mos, mock_10_multi_int_spec)

    # no int_times - only int_num is added to spec
    # It is set to the integration number for all spectra - no integrations in multislit data.
    assert_equal(mock_10_multi_int_spec.spec[0].spec_table["INT_NUM"], np.arange(1, 11))


def test_populate_time_keywords_multislit_table(
    mock_nirspec_mos, mock_nirspec_bots, mock_10_spec, log_watcher
):
    mock_nirspec_mos.meta.exposure.nints = 10
    mock_nirspec_mos.int_times = mock_nirspec_bots.int_times

    watcher = log_watcher("jwst.extract_1d.extract", message="Not using INT_TIMES table")
    ex.populate_time_keywords(mock_nirspec_mos, mock_10_spec)
    watcher.assert_seen()

    # int_times present but not used - no update
    assert "INT_NUM" not in mock_10_spec.spec[0].spec_table.columns.names


def test_populate_time_keywords_averaged(
    mock_nirspec_fs_one_slit, mock_nirspec_bots, mock_10_spec, log_watcher
):
    mock_nirspec_fs_one_slit.meta.exposure.nints = 10
    mock_nirspec_fs_one_slit.int_times = mock_nirspec_bots.int_times

    watcher = log_watcher("jwst.extract_1d.extract", message="Not using INT_TIMES table")
    ex.populate_time_keywords(mock_nirspec_fs_one_slit, mock_10_spec)
    watcher.assert_seen()

    # int_times not used - no update
    assert "INT_NUM" not in mock_10_spec.spec[0].spec_table.columns.names


def test_populate_time_keywords_mismatched_table(mock_nirspec_bots, mock_10_spec, log_watcher):
    # mock 20 integrations - table has 10
    mock_nirspec_bots.data = np.vstack([mock_nirspec_bots.data, mock_nirspec_bots.data])
    watcher = log_watcher("jwst.extract_1d.extract", message="Not using INT_TIMES table")
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_spec)
    watcher.assert_seen()

    # int_times not used - no update
    assert "INT_NUM" not in mock_10_spec.spec[0].spec_table.columns.names


def test_populate_time_keywords_missing_ints(mock_nirspec_bots, mock_10_spec, log_watcher):
    mock_nirspec_bots.meta.exposure.integration_start = 20
    watcher = log_watcher("jwst.extract_1d.extract", message="does not include rows")
    ex.populate_time_keywords(mock_nirspec_bots, mock_10_spec)
    watcher.assert_seen()

    # int_times not used - no update
    assert "INT_NUM" not in mock_10_spec.spec[0].spec_table.columns.names


def test_populate_time_keywords_ifu_table(
    mock_miri_ifu, mock_nirspec_bots, mock_10_spec, log_watcher
):
    mock_miri_ifu.meta.exposure.nints = 10
    mock_miri_ifu.int_times = mock_nirspec_bots.int_times

    watcher = log_watcher("jwst.extract_1d.extract", message="ignored for IFU")
    ex.populate_time_keywords(mock_miri_ifu, mock_10_spec)
    watcher.assert_seen()

    # int_times present but not used - no update
    assert "INT_NUM" not in mock_10_spec.spec[0].spec_table.columns.names


def test_populate_time_keywords_mismatched_spec(
    mock_nirspec_bots, mock_2_multi_int_spec, log_watcher
):
    watcher = log_watcher("jwst.extract_1d.extract", message="Don't understand n_output_spec")
    ex.populate_time_keywords(mock_nirspec_bots, mock_2_multi_int_spec)
    watcher.assert_seen()


def test_get_spectral_order(mock_nirspec_fs_one_slit):
    slit = mock_nirspec_fs_one_slit

    slit.meta.wcsinfo.spectral_order = 2
    assert ex.get_spectral_order(slit) == 2

    slit.meta.wcsinfo.spectral_order = None
    assert ex.get_spectral_order(slit) == 1

    slit.meta.wcsinfo = None
    assert ex.get_spectral_order(slit) == 1

    del slit.meta.wcsinfo
    assert ex.get_spectral_order(slit) == 1


def test_is_prism_false_conditions(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    assert ex.is_prism(model) is False

    model.meta.instrument.filter = None
    model.meta.instrument.grating = None
    assert ex.is_prism(model) is False

    model.meta.instrument.name = None
    assert ex.is_prism(model) is False


def test_is_prism_nirspec(mock_nirspec_fs_one_slit):
    mock_nirspec_fs_one_slit.meta.instrument.grating = "PRISM"
    assert ex.is_prism(mock_nirspec_fs_one_slit) is True


def test_is_prism_miri(mock_miri_ifu):
    mock_miri_ifu.meta.instrument.filter = "P750L"
    assert ex.is_prism(mock_miri_ifu) is True


def test_copy_keyword_info(mock_nirspec_fs_one_slit, mock_one_spec):
    expected = {
        "slitlet_id": 2,
        "source_id": 3,
        "source_name": "4",
        "source_alias": "5",
        "source_type": "POINT",
        "stellarity": 0.5,
        "source_xpos": -0.5,
        "source_ypos": -0.5,
        "source_ra": 10.0,
        "source_dec": 10.0,
        "shutter_state": "x",
        "wavelength_corrected": True,
        "pathloss_correction_type": "POINT",
        "barshadow_corrected": False,
    }
    for key, value in expected.items():
        setattr(mock_nirspec_fs_one_slit, key, value)
        assert not hasattr(mock_one_spec, key)

    ex.copy_keyword_info(mock_nirspec_fs_one_slit, "slit_name", mock_one_spec)
    assert mock_one_spec.name == "slit_name"
    for key, value in expected.items():
        assert getattr(mock_one_spec, key) == value


@pytest.mark.parametrize("partial", [True, False])
@pytest.mark.parametrize(
    "lower,upper",
    [
        (0, 19),
        (-1, 21),
        (np.full(10, 0.0), np.full(10, 19.0)),
        (np.linspace(-1, 0, 10), np.linspace(19, 20, 10)),
    ],
)
def test_set_weights_from_limits_whole_array(lower, upper, partial):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=float)
    yidx, _ = np.mgrid[: shape[0], : shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=partial)
    assert_equal(profile, 1.0)


@pytest.mark.parametrize(
    "lower,upper",
    [
        (10, 12),
        (9.5, 12.5),
        (np.linspace(9.5, 10, 10), np.linspace(12, 12.5, 10)),
    ],
)
def test_set_weights_from_limits_whole_pixel(lower, upper):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[: shape[0], : shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=False)
    assert_equal(profile[10:13], 1.0)


@pytest.mark.parametrize(
    "lower,upper",
    [
        (9.5, 12.5),
        (np.linspace(9.5, 10, 10), np.linspace(12, 12.5, 10)),
    ],
)
def test_set_weights_from_limits_partial_pixel(lower, upper):
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[: shape[0], : shape[1]]

    ex._set_weight_from_limits(profile, yidx, lower, upper, allow_partial=True)
    assert_allclose(profile[10:13], 1.0)
    assert_allclose(profile[9], 10 - lower)
    assert_allclose(profile[13], upper - 12)


def test_set_weights_from_limits_overlap():
    shape = (20, 10)
    profile = np.zeros(shape, dtype=np.float32)
    yidx, _ = np.mgrid[: shape[0], : shape[1]]

    # Set an aperture with partial pixel edges
    ex._set_weight_from_limits(profile, yidx, 9.5, 10.5, allow_partial=True)
    assert_allclose(profile[9], 0.5)
    assert_allclose(profile[11], 0.5)
    assert_allclose(profile[12], 0.0)

    # Set an overlapping region in the same profile
    ex._set_weight_from_limits(profile, yidx, 9.8, 11.5, allow_partial=True)

    # Higher weight from previous profile remains
    assert_allclose(profile[9], 0.5)

    # Previous partial pixel is now fully included
    assert_allclose(profile[11], 1.0)

    # New partial weight set on upper limit
    assert_allclose(profile[12], 0.5)


def test_box_profile_horizontal(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1

    # Exclude 2 pixels at left and right
    params["xstart"] = 1.5
    params["xstop"] = 7.5

    # Exclude 2 pixels at top and bottom, set half pixel weights
    # for another pixel at top and bottom edges
    params["ystart"] = 2.5
    params["ystop"] = 6.5
    profile = ex.box_profile(shape, extract_defaults, wl_array)

    # ystart/stop sets partial weights, xstart/stop sets whole pixels only
    assert_equal(profile[2:3, 3:8], 0.5)
    assert_equal(profile[7:8, 3:8], 0.5)
    assert_equal(profile[3:7, 3:8], 1.0)
    assert_equal(profile[:2], 0.0)
    assert_equal(profile[8:], 0.0)
    assert_equal(profile[:, :2], 0.0)
    assert_equal(profile[:, 8:], 0.0)


def test_box_profile_vertical(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 2

    # Exclude 2 pixels at "left" and "right" - in transposed aperture
    params["ystart"] = 1.5
    params["ystop"] = 7.5

    # Exclude 2 pixels at "top" and "bottom", set half pixel weights
    # for another pixel at top and bottom edges
    params["xstart"] = 2.5
    params["xstop"] = 6.5
    profile = ex.box_profile(shape, extract_defaults, wl_array)

    # xstart/stop sets partial weights, ystart/stop sets whole pixels only
    assert_equal(profile[3:8, 2:3], 0.5)
    assert_equal(profile[3:8, 7:8], 0.5)
    assert_equal(profile[3:8, 3:7], 1.0)

    assert_equal(profile[:2], 0.0)
    assert_equal(profile[8:], 0.0)
    assert_equal(profile[:, :2], 0.0)
    assert_equal(profile[:, 8:], 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_box_profile_bkg_coeff(extract_defaults, dispaxis):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = dispaxis

    # the definition for bkg_coeff is half a pixel off from start/stop definitions -
    # this will set the equivalent of start/stop 0-2, 7-9 -
    # 3 pixels at top and bottom of the array
    params["bkg_coeff"] = [[-0.5], [2.5], [6.5], [9.5]]

    profile, lower, upper = ex.box_profile(
        shape, extract_defaults, wl_array, coefficients="bkg_coeff", return_limits=True
    )
    if dispaxis == 2:
        profile = profile.T

    assert_equal(profile[:3], 1.0)
    assert_equal(profile[7:], 1.0)
    assert_equal(profile[3:7], 0.0)
    assert lower == 0.0
    assert upper == 9.0


def test_box_profile_bkg_coeff_median(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1
    params["bkg_fit"] = "median"

    # Attempt to set partial pixels at middle edges
    params["bkg_coeff"] = [[-0.5], [3.0], [6.0], [9.5]]

    profile, lower, upper = ex.box_profile(
        shape, extract_defaults, wl_array, coefficients="bkg_coeff", return_limits=True
    )

    # partial pixels are not allowed for fit type median - the profile is
    # set for whole pixels only
    assert_equal(profile[:3], 1.0)
    assert_equal(profile[7:], 1.0)
    assert_equal(profile[3:7], 0.0)
    assert lower == 0.0
    assert upper == 9.0


@pytest.mark.parametrize("swap_order", [False, True])
def test_box_profile_bkg_coeff_poly(extract_defaults, swap_order):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1
    params["bkg_fit"] = "poly"

    # Attempt to set partial pixels at middle edges
    if swap_order:
        # upper region first - this should make no difference.
        params["bkg_coeff"] = [[6.0], [9.5], [-0.5], [3.0]]
    else:
        params["bkg_coeff"] = [[-0.5], [3.0], [6.0], [9.5]]

    profile, lower, upper = ex.box_profile(
        shape, extract_defaults, wl_array, coefficients="bkg_coeff", return_limits=True
    )

    # partial pixels are allowed for fit type poly
    assert_equal(profile[:3], 1.0)
    assert_equal(profile[7:], 1.0)
    assert_equal(profile[3], 0.5)
    assert_equal(profile[6], 0.5)
    assert_equal(profile[4:6], 0.0)
    assert lower == 0.0
    assert upper == 9.0


@pytest.mark.parametrize("independent_var", ["pixel", "wavelength"])
def test_box_profile_src_coeff_constant(extract_defaults, independent_var):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1
    params["independent_var"] = independent_var

    # the definition for src_coeff is half a pixel off from start/stop definitions -
    # this will set the equivalent of start/stop 3-6, excluding
    # 3 pixels at top and bottom of the array
    params["src_coeff"] = [[2.5], [6.5]]

    profile, lower, upper = ex.box_profile(
        shape, extract_defaults, wl_array, coefficients="src_coeff", return_limits=True
    )
    assert_equal(profile[3:7], 1.0)
    assert_equal(profile[:3], 0.0)
    assert_equal(profile[7:], 0.0)
    assert lower == 3.0
    assert upper == 6.0


@pytest.mark.parametrize("independent_var", ["pixel", "wavelength"])
def test_box_profile_src_coeff_linear(extract_defaults, independent_var):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1
    params["independent_var"] = independent_var

    if independent_var == "wavelength":
        slope = 1 / (wl_array[0, 1] - wl_array[0, 0])
        start = -0.5 - wl_array[0, 0] * slope
        stop = start + 4
    else:
        slope = 1.0
        start = -0.5
        stop = 3.5

    # Set linearly increasing upper and lower edges,
    # starting at the bottom of the array, with a width of 4 pixels
    params["src_coeff"] = [[start, slope], [stop, slope]]
    expected = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
    ]

    profile, lower, upper = ex.box_profile(
        shape, extract_defaults, wl_array, coefficients="src_coeff", return_limits=True
    )
    assert_allclose(profile, expected, atol=1e-7)

    # upper and lower limits are averages
    assert_allclose(lower, 4.5)
    assert_allclose(upper, 7.5)


def test_box_profile_mismatched_coeff(extract_defaults):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = 1

    # set some mismatched coefficients limits
    params["src_coeff"] = [[2.5], [6.5], [9.5]]

    with pytest.raises(RuntimeError, match="must contain alternating lists"):
        ex.box_profile(shape, extract_defaults, wl_array)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_box_profile_from_width(extract_defaults, dispaxis):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = dispaxis

    if dispaxis == 1:
        # Set ystart and ystop to center on pixel 4
        params["ystart"] = 2.0
        params["ystop"] = 6.0
    else:
        # Set xstart and xstop to center on pixel 4
        params["xstart"] = 2.0
        params["xstop"] = 6.0

    # Set a width to 6 pixels
    params["extract_width"] = 6.0

    profile = ex.box_profile(shape, extract_defaults, wl_array)
    if dispaxis == 2:
        profile = profile.T

    # Aperture is centered at pixel 4, to start at 1.5, end at 6.5
    assert_equal(profile[2:7], 1.0)
    assert_equal(profile[1], 0.5)
    assert_equal(profile[7], 0.5)
    assert_equal(profile[0], 0.0)
    assert_equal(profile[8:], 0.0)


@pytest.mark.parametrize("dispaxis", [1, 2])
def test_box_profile_from_trace(extract_defaults, dispaxis):
    shape = (10, 10)
    wl_array = np.empty(shape)
    wl_array[:] = np.linspace(3, 5, 10)

    params = extract_defaults
    params["dispaxis"] = dispaxis

    # Set a linear trace
    params["trace"] = np.arange(10) + 1.5

    # Set the width to 4 pixels
    params["extract_width"] = 4.0

    # Make the profile
    profile, lower, upper = ex.box_profile(shape, extract_defaults, wl_array, return_limits=True)
    if dispaxis == 2:
        profile = profile.T

    expected = [
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
    ]

    assert_allclose(profile, expected)

    # upper and lower limits are averages
    assert_allclose(lower, 4.5)
    assert_allclose(upper, 7.5)


@pytest.mark.parametrize("middle", [None, 7])
@pytest.mark.parametrize("dispaxis", [1, 2])
def test_aperture_center(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[1:4] = 1.0
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 2.0
    if middle is None:
        assert spec_center == 4.5
    else:
        assert spec_center == middle


@pytest.mark.parametrize("middle", [None, 7])
@pytest.mark.parametrize("dispaxis", [1, 2])
def test_aperture_center_zero_weight(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    slit_center, spec_center = ex.aperture_center(profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 4.5
    if middle is None:
        assert spec_center == 4.5
    else:
        assert spec_center == middle


@pytest.mark.parametrize("middle", [None, 7])
@pytest.mark.parametrize("dispaxis", [1, 2])
def test_aperture_center_variable_weight_by_slit(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[1:4] = np.arange(10)
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(profile, dispaxis=dispaxis, middle_pix=middle)
    assert slit_center == 2.0
    if middle is None:
        assert_allclose(spec_center, 6.3333333)
    else:
        assert spec_center == middle


@pytest.mark.parametrize("middle", [None, 7])
@pytest.mark.parametrize("dispaxis", [1, 2])
def test_aperture_center_variable_weight_by_spec(middle, dispaxis):
    profile = np.zeros((10, 10), dtype=np.float32)
    profile[:, 1:4] = np.arange(10)[:, None]
    if dispaxis != 1:
        profile = profile.T
    slit_center, spec_center = ex.aperture_center(profile, dispaxis=dispaxis, middle_pix=middle)
    if middle is None:
        assert_allclose(slit_center, 6.3333333)
        assert_allclose(spec_center, 2.0)
    else:
        assert_allclose(slit_center, 4.5)
        assert spec_center == middle


def test_shift_by_offset_horizontal(extract_defaults):
    offset = 2.5

    extract_params = extract_defaults.copy()
    extract_params["dispaxis"] = 1
    extract_params["position_offset"] = offset

    ex.shift_by_offset(offset, extract_params)
    assert extract_params["xstart"] == extract_defaults["xstart"]
    assert extract_params["xstop"] == extract_defaults["xstop"]
    assert extract_params["ystart"] == extract_defaults["ystart"] + offset
    assert extract_params["ystop"] == extract_defaults["ystop"] + offset


def test_shift_by_offset_vertical(extract_defaults):
    offset = 2.5

    extract_params = extract_defaults.copy()
    extract_params["dispaxis"] = 2
    extract_params["position_offset"] = offset

    ex.shift_by_offset(offset, extract_params)
    assert extract_params["xstart"] == extract_defaults["xstart"] + offset
    assert extract_params["xstop"] == extract_defaults["xstop"] + offset
    assert extract_params["ystart"] == extract_defaults["ystart"]
    assert extract_params["ystop"] == extract_defaults["ystop"]


def test_shift_by_offset_coeff(extract_defaults):
    offset = 2.5

    extract_params = extract_defaults.copy()
    extract_params["dispaxis"] = 1
    extract_params["position_offset"] = offset
    extract_params["src_coeff"] = [[2.5, 1.0], [6.5, 1.0]]
    extract_params["bkg_coeff"] = [[-0.5], [3.0], [6.0], [9.5]]

    ex.shift_by_offset(offset, extract_params)
    assert extract_params["src_coeff"] == [[2.5 + offset, 1.0], [6.5 + offset, 1.0]]
    assert extract_params["bkg_coeff"] == [
        [-0.5 + offset],
        [3.0 + offset],
        [6.0 + offset],
        [9.5 + offset],
    ]


def test_shift_by_offset_trace(extract_defaults):
    offset = 2.5

    extract_params = extract_defaults.copy()
    extract_params["dispaxis"] = 1
    extract_params["position_offset"] = offset
    extract_params["trace"] = np.arange(10, dtype=float)

    ex.shift_by_offset(offset, extract_params, update_trace=True)
    assert_equal(extract_params["trace"], np.arange(10) + offset)


def test_shift_by_offset_trace_no_update(extract_defaults):
    offset = 2.5

    extract_params = extract_defaults.copy()
    extract_params["dispaxis"] = 1
    extract_params["position_offset"] = offset
    extract_params["trace"] = np.arange(10, dtype=float)

    ex.shift_by_offset(offset, extract_params, update_trace=False)
    assert_equal(extract_params["trace"], np.arange(10))


@pytest.mark.parametrize("is_slit", [True, False])
def test_define_aperture_nirspec(mock_nirspec_fs_one_slit, extract_defaults, is_slit):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    if is_slit:
        slit = model
    else:
        slit = None
    exptype = "NRS_FIXEDSLIT"
    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec, wavelength, profile, bg_profile, nod_profile, limits = result
    assert_allclose(ra, 45.05)
    assert_allclose(dec, 45.1)
    assert wavelength.shape == (model.data.shape[1],)
    assert profile.shape == model.data.shape

    # Default profile is the full array
    assert_equal(profile, 1.0)
    assert limits == (0, model.data.shape[0] - 1, 0, model.data.shape[1] - 1)

    # Default bg profile is None
    assert bg_profile is None


@pytest.mark.parametrize("is_slit", [True, False])
def test_define_aperture_miri(mock_miri_lrs_fs, extract_defaults, is_slit):
    model = mock_miri_lrs_fs
    extract_defaults["dispaxis"] = 2
    if is_slit:
        slit = model
    else:
        slit = None
    exptype = "MIR_LRS-FIXEDSLIT"
    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec, wavelength, profile, bg_profile, nod_profile, limits = result
    assert_allclose(ra, 45.05)
    assert_allclose(dec, 45.1)
    assert wavelength.shape == (model.data.shape[1],)
    assert profile.shape == model.data.shape

    # Default profile is the full array
    assert_equal(profile, 1.0)
    assert limits == (0, model.data.shape[0] - 1, 0, model.data.shape[1] - 1)

    # Default bg profile is None
    assert bg_profile is None


def test_define_aperture_with_bg(mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    slit = None
    exptype = "NRS_FIXEDSLIT"

    extract_defaults["subtract_background"] = True
    extract_defaults["bkg_coeff"] = [[-0.5], [2.5]]

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    bg_profile = result[-3]

    # Bg profile has 1s in the first 3 rows
    assert bg_profile.shape == model.data.shape
    assert_equal(bg_profile[:3], 1.0)
    assert_equal(bg_profile[3:], 0.0)


def test_define_aperture_empty_aperture(mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    slit = None
    exptype = "NRS_FIXEDSLIT"

    # Set the extraction limits out of range
    extract_defaults["ystart"] = 2000
    extract_defaults["ystop"] = 3000

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, _, _, limits = result

    assert_equal(profile, 0.0)
    assert limits == (2000, 3000, None, None)


def test_define_aperture_bad_wcs(monkeypatch, mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    slit = None
    exptype = "NRS_FIXEDSLIT"

    # Set a wavelength so wcs is not called to retrieve it
    model.wavelength = np.empty_like(model.data)
    model.wavelength[:] = np.linspace(3, 5, model.data.shape[1])

    # mock a bad wcs
    def return_nan(*args):
        return np.nan, np.nan, np.nan

    monkeypatch.setattr(model.meta, "wcs", return_nan)

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    ra, dec = result[:2]

    # RA and Dec returned are none
    assert ra is None
    assert dec is None


def test_define_aperture_use_source(monkeypatch, mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    slit = None
    exptype = "NRS_FIXEDSLIT"

    # mock the source location function
    def mock_source_location(*args):
        return 24, 7.74, 9.5, np.full(model.data.shape[-1], 9.5)

    monkeypatch.setattr(ex, "location_from_wcs", mock_source_location)

    # set parameters to extract a 6 pixel aperture, centered on source location
    extract_defaults["use_source_posn"] = True
    extract_defaults["extract_width"] = 6.0

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, _, _, limits = result

    assert_equal(profile[:7], 0.0)
    assert_equal(profile[7:13], 1.0)
    assert_equal(profile[13:], 0.0)


def test_define_aperture_extra_offset(mock_nirspec_fs_one_slit, extract_defaults):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    slit = None
    exptype = "NRS_FIXEDSLIT"

    extract_defaults["position_offset"] = 2.0

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, _, _, limits = result
    assert profile.shape == model.data.shape

    # Default profile is shifted 2 pixels up
    assert_equal(profile[:2], 0.0)
    assert_equal(profile[2:], 1.0)
    assert limits == (2, model.data.shape[0] + 1, 0, model.data.shape[1] - 1)


def test_define_aperture_optimal(mock_miri_lrs_fs, extract_defaults, psf_reference_file):
    model = mock_miri_lrs_fs
    extract_defaults["dispaxis"] = 2
    slit = None
    exptype = "MIR_LRS-FIXEDSLIT"

    # set parameters for optimal extraction
    extract_defaults["extraction_type"] = "optimal"
    extract_defaults["use_source_posn"] = True
    extract_defaults["psf"] = psf_reference_file

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, bg_profile, nod_profile, limits = result

    assert bg_profile is None
    assert nod_profile is None

    # profile is normalized along cross-dispersion
    assert_allclose(np.sum(profile, axis=1), 1.0)

    # trace is centered on 1.0, near the edge of the slit,
    # and the psf data has the same size as the array (50x50),
    # so only half the psf is included
    npix = 26
    assert_equal(np.sum(profile != 0, axis=1), npix)

    # psf is uniform when in range, 0 otherwise
    assert_allclose(profile[:, :npix], 1 / npix)
    assert_allclose(profile[:, npix:], 0.0)


def test_define_aperture_optimal_with_nod(
    monkeypatch, mock_miri_lrs_fs, extract_defaults, psf_reference_file
):
    model = mock_miri_lrs_fs
    extract_defaults["dispaxis"] = 2
    slit = None
    exptype = "MIR_LRS-FIXEDSLIT"

    # mock nod subtraction
    mock_miri_lrs_fs.meta.cal_step.bkg_subtract = "COMPLETE"
    mock_miri_lrs_fs.meta.dither.primary_type = "ALONG-SLIT-NOD"

    # mock a nod position at the opposite end of the array
    def mock_nod(*args, **kwargs):
        return 48.0

    monkeypatch.setattr(pp, "nod_pair_location", mock_nod)

    # set parameters for optimal extraction
    extract_defaults["extraction_type"] = "optimal"
    extract_defaults["use_source_posn"] = True
    extract_defaults["psf"] = psf_reference_file
    extract_defaults["model_nod_pair"] = True

    result = ex.define_aperture(model, slit, extract_defaults, exptype)
    _, _, _, profile, bg_profile, nod_profile, limits = result

    assert bg_profile is None
    assert nod_profile is not None

    # profiles are normalized along cross-dispersion,
    # nod profile is negative
    assert_allclose(np.sum(profile, axis=1), 1.0)
    assert_allclose(np.sum(nod_profile, axis=1), -1.0)

    # positive trace is centered on 1.0, negative trace on
    # 48.0, array size is 50.
    npix = 26
    assert_equal(np.sum(profile != 0, axis=1), npix)
    assert_equal(np.sum(nod_profile != 0, axis=1), npix)

    # psf is uniform when in range, 0 otherwise
    assert_allclose(profile[:, :npix], 1 / npix)
    assert_allclose(profile[:, npix:], 0.0)

    # nod profile is the same, but negative, and at the other
    # end of the array
    assert_allclose(nod_profile[:, -npix:], -1 / npix)
    assert_allclose(nod_profile[:, :-npix], 0.0)


def test_extract_one_slit_horizontal(
    mock_nirspec_fs_one_slit, extract_defaults, simple_profile, background_profile
):
    # update parameters to subtract background
    extract_defaults["dispaxis"] = 1
    extract_defaults["subtract_background"] = True
    extract_defaults["bkg_fit"] = "poly"
    extract_defaults["bkg_order"] = 1

    # set a source in the profile region
    mock_nirspec_fs_one_slit.data[simple_profile != 0] += 1.0

    result = ex.extract_one_slit(
        mock_nirspec_fs_one_slit, -1, simple_profile, background_profile, None, extract_defaults
    )

    for data in result[:-2]:
        assert np.all(data > 0)
        assert data.shape == (mock_nirspec_fs_one_slit.data.shape[1],)

    # residuals from the 2D scene model should be zero - this simple case
    # is exactly modeled with a box profile
    scene_model = result[-2]
    assert scene_model.shape == mock_nirspec_fs_one_slit.data.shape
    assert_allclose(np.abs(mock_nirspec_fs_one_slit.data - scene_model), 0, atol=1e-7)

    residual = result[-1]
    assert residual.shape == mock_nirspec_fs_one_slit.data.shape
    assert_allclose(np.abs(residual), 0, atol=1e-7)

    # flux should be 1.0 * npixels
    flux = result[0]
    npixels = result[-3]
    assert_allclose(flux, npixels)

    # npixels is sum of profile
    assert_equal(npixels, np.sum(simple_profile, axis=0))


def test_extract_one_slit_vertical(
    mock_miri_lrs_fs, extract_defaults, simple_profile, background_profile
):
    model = mock_miri_lrs_fs
    profile = simple_profile.T
    profile_bg = background_profile.T

    # update parameters to subtract background
    extract_defaults["dispaxis"] = 2
    extract_defaults["subtract_background"] = True
    extract_defaults["bkg_fit"] = "poly"
    extract_defaults["bkg_order"] = 1

    # set a source in the profile region
    model.data[profile != 0] += 1.0

    result = ex.extract_one_slit(model, -1, profile, profile_bg, None, extract_defaults)

    for data in result[:-2]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # residuals from the 2D scene model should be zero - this simple case
    # is exactly modeled with a box profile
    scene_model = result[-2]
    assert scene_model.shape == model.data.shape
    assert_allclose(np.abs(model.data - scene_model), 0, atol=1e-7)

    residual = result[-1]
    assert residual.shape == model.data.shape
    assert_allclose(np.abs(residual), 0, atol=1e-7)

    # flux should be 1.0 * npixels
    flux = result[0]
    npixels = result[-3]
    assert_allclose(flux, npixels)

    # npixels is sum of profile
    assert_equal(npixels, np.sum(profile, axis=1))


def test_extract_one_slit_vertical_no_bg(mock_miri_lrs_fs, extract_defaults, simple_profile):
    model = mock_miri_lrs_fs
    profile = simple_profile.T
    extract_defaults["dispaxis"] = 2

    result = ex.extract_one_slit(model, -1, profile, None, None, extract_defaults)

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[0],)

    # npixels is the sum of the profile
    assert_allclose(result[8], np.sum(simple_profile, axis=0))

    # scene model and residual has 2D shape
    assert result[-2].shape == model.data.shape
    assert result[-1].shape == model.data.shape


def test_extract_one_slit_multi_int(
    mock_nirspec_bots, extract_defaults, simple_profile, log_watcher
):
    model = mock_nirspec_bots
    extract_defaults["dispaxis"] = 1

    watcher = log_watcher("jwst.extract_1d.extract", message="Extracting integration 2")
    result = ex.extract_one_slit(model, 1, simple_profile, None, None, extract_defaults)
    watcher.assert_seen()

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[2],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[2],)

    # npixels is the sum of the profile
    assert_allclose(result[8], np.sum(simple_profile, axis=0))

    # scene model and residual has 2D shape
    assert result[-2].shape == model.data.shape[-2:]
    assert result[-1].shape == model.data.shape[-2:]


def test_extract_one_slit_missing_var(mock_nirspec_fs_one_slit, extract_defaults, simple_profile):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1

    # Test that mismatched variances still extract okay.
    # This is probably only possible for var_flat, which is optional and
    # uninitialized if flat fielding is skipped, but the code has handling
    # for all 3 variance arrays.
    model.var_rnoise = np.zeros((10, 10))
    model.var_poisson = np.zeros((10, 10))
    model.var_flat = np.zeros((10, 10))

    result = ex.extract_one_slit(model, -1, simple_profile, None, None, extract_defaults)

    # flux is nonzero
    assert np.all(result[0] > 0)
    assert result[0].shape == (model.data.shape[1],)

    # variances are zero
    for data in result[1:4]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[1],)


def test_extract_one_slit_optimal_horizontal(
    mock_nirspec_fs_one_slit, extract_defaults, nod_profile, negative_nod_profile
):
    model = mock_nirspec_fs_one_slit
    extract_defaults["dispaxis"] = 1
    extract_defaults["extraction_type"] = "optimal"

    result = ex.extract_one_slit(
        model, -1, nod_profile, None, negative_nod_profile, extract_defaults
    )

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[0],)

    # npixels is the sum of the pixels in the positive profile
    assert_allclose(result[8], np.sum(nod_profile > 0, axis=0))

    # scene model and residual has 2D shape
    assert result[-1].shape == model.data.shape
    assert result[-2].shape == model.data.shape


def test_extract_one_slit_optimal_vertical(
    mock_miri_lrs_fs, extract_defaults, nod_profile, negative_nod_profile
):
    model = mock_miri_lrs_fs
    nod_profile = nod_profile.T
    negative_nod_profile = negative_nod_profile.T
    extract_defaults["dispaxis"] = 2
    extract_defaults["extraction_type"] = "optimal"

    result = ex.extract_one_slit(
        model, -1, nod_profile, None, negative_nod_profile, extract_defaults
    )

    # flux and variances are nonzero
    for data in result[:4]:
        assert np.all(data > 0)
        assert data.shape == (model.data.shape[0],)

    # background and variances are zero
    for data in result[4:8]:
        assert np.all(data == 0)
        assert data.shape == (model.data.shape[0],)

    # npixels is the sum of the pixels in the positive profile
    assert_allclose(result[8], np.sum(nod_profile > 0, axis=1))

    # scene model and residual has 2D shape
    assert result[-2].shape == model.data.shape
    assert result[-1].shape == model.data.shape


def test_create_extraction_with_photom(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.meta.cal_step.photom = "COMPLETE"

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].spec_table.columns["flux"].unit == "Jy"


def test_create_extraction_without_photom(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.meta.cal_step.photom = "SKIPPED"

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].spec_table.columns["flux"].unit == "DN/s"


def test_create_extraction_missing_src_type(create_extraction_inputs):
    model = create_extraction_inputs[0]
    model.source_type = None
    model.meta.target.source_type = "EXTENDED"

    ex.create_extraction(*create_extraction_inputs)

    output_model = create_extraction_inputs[2]
    assert output_model.spec[0].source_type == "EXTENDED"


def test_create_extraction_no_match(create_extraction_inputs):
    create_extraction_inputs[4] = "bad slitname"
    with pytest.raises(ValueError, match="Missing extraction parameters"):
        ex.create_extraction(*create_extraction_inputs)


def test_create_extraction_partial_match(create_extraction_inputs, log_watcher):
    # match a slit that has a mismatched spectral order specified
    create_extraction_inputs[4] = "slit7"

    watcher = log_watcher("jwst.extract_1d.extract", message="Spectral order 1 not found")
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    watcher.assert_seen()


def test_create_extraction_missing_dispaxis(create_extraction_inputs, log_watcher):
    create_extraction_inputs[0].meta.wcsinfo.dispersion_direction = None
    watcher = log_watcher(
        "jwst.extract_1d.extract", message="dispersion direction information is missing"
    )
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    watcher.assert_seen()


def test_create_extraction_missing_wavelengths(create_extraction_inputs, log_watcher):
    model = create_extraction_inputs[0]
    model.wavelength = np.full_like(model.data, np.nan)
    watcher = log_watcher("jwst.extract_1d.extract", message="Spectrum is empty; no valid data")
    with pytest.raises(ex.ContinueError):
        ex.create_extraction(*create_extraction_inputs)
    watcher.assert_seen()


def test_create_extraction_nrs_apcorr(
    create_extraction_inputs, nirspec_fs_apcorr, mock_nirspec_bots, log_watcher
):
    model = mock_nirspec_bots
    model.source_type = "POINT"
    model.meta.cal_step.photom = "COMPLETE"
    create_extraction_inputs[0] = model

    watcher = log_watcher("jwst.extract_1d.extract", message="Tabulating aperture correction")
    ex.create_extraction(
        *create_extraction_inputs, apcorr_ref_model=nirspec_fs_apcorr, use_source_posn=False
    )
    watcher.assert_seen()


def test_create_extraction_one_int(create_extraction_inputs, mock_nirspec_bots, log_watcher):
    # Input model is a cube, but with only one integration
    model = mock_nirspec_bots
    model.data = model.data[0].reshape(1, *model.data.shape[-2:])
    create_extraction_inputs[0] = model
    create_extraction_inputs[4] = "S1600A1"

    watcher = log_watcher("jwst.extract_1d.extract", message="1 integration done")
    ex.create_extraction(*create_extraction_inputs, log_increment=1)
    output_model = create_extraction_inputs[2]
    assert len(output_model.spec) == 1
    watcher.assert_seen()


def test_create_extraction_log_increment(create_extraction_inputs, mock_nirspec_bots, log_watcher):
    create_extraction_inputs[0] = mock_nirspec_bots
    create_extraction_inputs[4] = "S1600A1"

    # all integrations are logged
    watcher = log_watcher("jwst.extract_1d.extract", message="... 9 integrations done")
    ex.create_extraction(*create_extraction_inputs, log_increment=1)
    watcher.assert_seen()


@pytest.mark.parametrize("use_source", [True, False, None])
@pytest.mark.parametrize("source_type", ["POINT", "EXTENDED"])
def test_create_extraction_use_source(
    monkeypatch,
    create_extraction_inputs,
    mock_nirspec_fs_one_slit,
    use_source,
    source_type,
    log_watcher,
):
    model = mock_nirspec_fs_one_slit
    model.source_type = source_type
    create_extraction_inputs[0] = model

    # mock the source location function
    def mock_source_location(*args):
        return 24, 7.74, 9.5, np.full(model.data.shape[-1], 9.5)

    monkeypatch.setattr(ex, "location_from_wcs", mock_source_location)

    watcher = log_watcher("jwst.extract_1d.extract")
    if source_type != "POINT" and use_source is None:
        # If not specified, source position should be used only if POINT
        watcher.message = "Setting use_source_posn to False"
    elif use_source is True or (source_type == "POINT" and use_source is None):
        # If explicitly set to True, or unspecified + source type is POINT,
        # source position is used
        watcher.message = "Aperture start/stop: -15"
    else:
        # If False, source position is not used
        watcher.message = "Aperture start/stop: 0"
    ex.create_extraction(*create_extraction_inputs, use_source_posn=use_source)
    watcher.assert_seen()


@pytest.mark.parametrize("extract_width", [None, 7])
def test_create_extraction_use_trace(
    monkeypatch, create_extraction_inputs, mock_nirspec_bots, extract_width, log_watcher
):
    model = mock_nirspec_bots
    create_extraction_inputs[0] = model
    aper = create_extraction_inputs[3]["apertures"]
    create_extraction_inputs[4] = "S1600A1"
    for i in range(len(aper)):
        if aper[i]["id"] == "S1600A1":
            aper[i]["use_source_posn"] = True
            aper[i]["extract_width"] = extract_width
            aper[i]["position_offset"] = 0

    # mock the source location function
    def mock_source_location(*args):
        return 24, 7.74, 25, np.full(model.data.shape[-1], 25)

    monkeypatch.setattr(ex, "location_from_wcs", mock_source_location)
    watcher = log_watcher("jwst.extract_1d.extract")
    if extract_width is not None:
        # If explicitly set to True, or unspecified + source type is POINT,
        # source position is used
        watcher.message = "aperture start/stop from trace: 22"
    else:
        # If False, source trace is not used
        watcher.message = "Aperture start/stop: 0"
    ex.create_extraction(*create_extraction_inputs)
    watcher.assert_seen()


def test_create_extraction_optimal(monkeypatch, create_extraction_inputs, psf_reference_file):
    model = create_extraction_inputs[0]

    # mock nod subtraction
    model.meta.cal_step.bkg_subtract = "COMPLETE"
    model.meta.dither.primary_type = "2-POINT-NOD"

    # mock a nod position at the opposite end of the array
    def mock_nod(*args, **kwargs):
        return 48.0

    monkeypatch.setattr(pp, "nod_pair_location", mock_nod)

    profile_model, _, _ = ex.create_extraction(
        *create_extraction_inputs,
        save_profile=True,
        psf_ref_name=psf_reference_file,
        use_source_posn=True,
        extraction_type="optimal",
        model_nod_pair=True,
    )

    assert profile_model is not None

    # profile contains positive and negative nod, summed
    assert np.all(profile_model.data[:10] > 0)
    assert np.all(profile_model.data[-10:] < 0)

    profile_model.close()


def test_run_extract1d(mock_nirspec_mos):
    model = mock_nirspec_mos
    output_model, profile_model, scene_model, residual = ex.run_extract1d(model)
    assert isinstance(output_model, dm.MultiSpecModel)
    assert profile_model is None
    assert scene_model is None
    assert residual is None
    output_model.close()


def test_run_extract1d_save_models(mock_niriss_wfss_l3):
    model = mock_niriss_wfss_l3
    output_model, profile_model, scene_model, residual = ex.run_extract1d(
        model, save_profile=True, save_scene_model=True, save_residual_image=True
    )
    assert isinstance(output_model, dm.MultiSpecModel)
    assert isinstance(profile_model, ModelContainer)
    assert isinstance(scene_model, ModelContainer)
    assert isinstance(residual, ModelContainer)

    assert len(profile_model) == len(model)
    assert len(scene_model) == len(model)
    assert len(residual) == len(model)

    for pmodel in profile_model:
        assert isinstance(pmodel, dm.ImageModel)
    for smodel in scene_model:
        assert isinstance(smodel, dm.ImageModel)
    for rmodel in residual:
        assert isinstance(rmodel, dm.ImageModel)

    output_model.close()
    profile_model.close()
    scene_model.close()
    residual.close()


def test_run_extract1d_save_cube_scene(mock_nirspec_bots):
    model = mock_nirspec_bots
    output_model, profile_model, scene_model, residual = ex.run_extract1d(
        model, save_profile=True, save_scene_model=True, save_residual_image=True
    )
    assert isinstance(output_model, dm.TSOMultiSpecModel)
    assert isinstance(profile_model, dm.ImageModel)
    assert isinstance(scene_model, dm.CubeModel)
    assert isinstance(residual, dm.CubeModel)

    assert profile_model.data.shape == model.data.shape[-2:]
    assert scene_model.data.shape == model.data.shape
    assert residual.data.shape == model.data.shape

    output_model.close()
    profile_model.close()
    scene_model.close()
    residual.close()


def test_run_extract1d_tso(mock_nirspec_bots):
    model = mock_nirspec_bots
    output_model, _, _, _ = ex.run_extract1d(model)

    # time and integration keywords are populated
    for i, spec in enumerate(output_model.spec[0].spec_table):
        assert spec["int_num"] == i + 1
        assert_allclose(spec["MJD-BEG"], 59729.04367729)

    output_model.close()


@pytest.mark.parametrize("from_name_attr", [True, False])
def test_run_extract1d_slitmodel_name(mock_nirspec_fs_one_slit, from_name_attr):
    model = mock_nirspec_fs_one_slit
    slit_name = "S200A1"
    if from_name_attr:
        model.name = slit_name
        model.meta.instrument.fixed_slit = None
    else:
        model.name = None
        model.meta.instrument.fixed_slit = "S200A1"

    output_model, _, _, _ = ex.run_extract1d(model)
    assert output_model.spec[0].name == "S200A1"

    output_model.close()


@pytest.mark.parametrize("from_name_attr", [True, False])
def test_run_extract1d_imagemodel_name(mock_miri_lrs_fs, from_name_attr):
    model = mock_miri_lrs_fs
    slit_name = "test_slit_name"
    if from_name_attr:
        model.name = slit_name
    else:
        model.name = None

    output_model, _, _, _ = ex.run_extract1d(model)
    if from_name_attr:
        assert output_model.spec[0].name == "test_slit_name"
    else:
        assert output_model.spec[0].name == "MIR_LRS-FIXEDSLIT"
    output_model.close()


def test_run_extract1d_apcorr(mock_miri_lrs_fs, miri_lrs_apcorr_file, log_watcher):
    model = mock_miri_lrs_fs
    model.meta.target.source_type = "POINT"

    watcher = log_watcher("jwst.extract_1d.extract", message="Creating aperture correction")
    output_model, _, _, _ = ex.run_extract1d(model, apcorr_ref_name=miri_lrs_apcorr_file)
    watcher.assert_seen()

    output_model.close()


def test_run_extract1d_apcorr_optimal(
    mock_miri_lrs_fs, miri_lrs_apcorr_file, psf_reference_file, log_watcher
):
    model = mock_miri_lrs_fs
    model.meta.target.source_type = "POINT"

    # Aperture correction that is otherwise valid is nonetheless
    # turned off for optimal extraction
    watcher = log_watcher(
        "jwst.extract_1d.extract", message="Turning off aperture correction for optimal extraction"
    )
    output_model, _, _, _ = ex.run_extract1d(
        model,
        apcorr_ref_name=miri_lrs_apcorr_file,
        psf_ref_name=psf_reference_file,
        extraction_type="optimal",
    )
    watcher.assert_seen()

    output_model.close()


def test_run_extract1d_optimal_no_psf(mock_miri_lrs_fs, log_watcher):
    model = mock_miri_lrs_fs
    model.meta.target.source_type = "POINT"

    # Optimal extraction is turned off if there is no psf file provided
    watcher = log_watcher("jwst.extract_1d.extract", message="Optimal extraction is not available")
    output_model, _, _, _ = ex.run_extract1d(model, extraction_type="optimal")
    watcher.assert_seen()

    output_model.close()


def test_run_extract1d_invalid():
    model = dm.MultiSpecModel()
    with pytest.raises(TypeError, match="Can't extract a spectrum"):
        ex.run_extract1d(model)


def test_run_extract1d_zeroth_order_slit(mock_nirspec_fs_one_slit):
    model = mock_nirspec_fs_one_slit
    model.meta.wcsinfo.spectral_order = 0
    output_model, _, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_zeroth_order_image(mock_miri_lrs_fs):
    model = mock_miri_lrs_fs
    model.meta.wcsinfo.spectral_order = 0
    model.meta.instrument.filter = None
    output_model, _, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_zeroth_order_multispec(mock_nirspec_mos):
    model = mock_nirspec_mos
    for slit in model.slits:
        slit.meta.wcsinfo.spectral_order = 0
    output_model, _, _, _ = ex.run_extract1d(model)

    # no spectra extracted for zeroth order
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_no_data(mock_niriss_wfss_l3):
    container = mock_niriss_wfss_l3
    for model in container:
        model.data = np.array([])
    output_model, _, _, _ = ex.run_extract1d(container)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_slit(monkeypatch, mock_nirspec_fs_one_slit):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError("Test error")

    monkeypatch.setattr(ex, "create_extraction", raise_continue_error)
    output_model, _, _, _ = ex.run_extract1d(mock_nirspec_fs_one_slit)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_image(monkeypatch, mock_miri_lrs_fs):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError("Test error")

    monkeypatch.setattr(ex, "create_extraction", raise_continue_error)
    output_model, _, _, _ = ex.run_extract1d(mock_miri_lrs_fs)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()


def test_run_extract1d_continue_error_multislit(monkeypatch, mock_nirspec_mos):
    def raise_continue_error(*args, **kwargs):
        raise ex.ContinueError("Test error")

    monkeypatch.setattr(ex, "create_extraction", raise_continue_error)
    output_model, _, _, _ = ex.run_extract1d(mock_nirspec_mos)

    # no spectra extracted
    assert len(output_model.spec) == 0
    output_model.close()
