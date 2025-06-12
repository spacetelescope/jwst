import json
import pytest
import numpy as np
from jwst.badpix_selfcal.badpix_selfcal_step import BadpixSelfcalStep, _parse_inputs
from jwst.badpix_selfcal.badpix_selfcal import badpix_selfcal, apply_flags
from jwst import datamodels as dm
from jwst.assign_wcs import AssignWcsStep, miri
from gwcs import wcs
from astropy.io import fits
from stdatamodels.jwst.datamodels.dqflags import pixel

wcs_kw = {
    "wcsaxes": 3,
    "ra_ref": 165,
    "dec_ref": 54,
    "v2_ref": -8.3942412,
    "v3_ref": -5.3123744,
    "roll_ref": 37,
    "crpix1": 1024,
    "crpix2": 1024,
    "crpix3": 0,
    "cdelt1": 0.08,
    "cdelt2": 0.08,
    "cdelt3": 1,
    "ctype1": "RA---TAN",
    "ctype2": "DEC--TAN",
    "ctype3": "WAVE",
    "pc1_1": 1,
    "pc1_2": 0,
    "pc1_3": 0,
    "pc2_1": 0,
    "pc2_2": 1,
    "pc2_3": 0,
    "pc3_1": 0,
    "pc3_2": 0,
    "pc3_3": 1,
    "cunit1": "deg",
    "cunit2": "deg",
    "cunit3": "um",
}

# these are often warm pixels, so need to be much larger than noise,
# which has std dev of 1,
# but smaller than smoothly-varying background, which maxes around 20
hotpixel_intensity = 10
outlier_indices = [(100, 100), (300, 300), (500, 600), (1000, 900)]


def create_hdul(detector, channel, band):
    hdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdu.header["telescop"] = "JWST"
    phdu.header["filename"] = "test" + channel + band
    phdu.header["instrume"] = "MIRI"
    phdu.header["detector"] = detector
    phdu.header["CHANNEL"] = channel
    phdu.header["BAND"] = band
    phdu.header["time-obs"] = "8:59:37"
    phdu.header["date-obs"] = "2017-09-05"
    phdu.header["exp_type"] = "MIR_MRS"
    scihdu = fits.ImageHDU()
    scihdu.header["EXTNAME"] = "SCI"
    scihdu.header.update(wcs_kw)
    hdul.append(phdu)
    hdul.append(scihdu)
    return hdul


def create_reference_files(datamodel):
    refs = {}
    step = AssignWcsStep()
    for reftype in AssignWcsStep.reference_file_types:
        refs[reftype] = step.get_reference_file(datamodel, reftype)

    return refs


@pytest.fixture(scope="module")
def background():
    """
    Create background IFUImageModel for testing. This is a mockup of the expected
    background data in a .rate file.
    Three components: random noise, low-order variability, and outliers
    """

    # random noise
    rng = np.random.default_rng(seed=77)
    shp = (1024, 1032)
    noise = rng.standard_normal(shp)

    # make 2-d polynomial representing background level
    c = np.array([[1, 3, 5], [2, 4, 6]])
    x = np.linspace(-1, 1, shp[0])
    y = np.linspace(-1, 1, shp[1])
    low_order_variability = np.polynomial.polynomial.polygrid2d(x, y, c)

    # add some outliers
    outliers = np.zeros(shp)
    for idx in outlier_indices:
        outliers[idx] = hotpixel_intensity
    # one negative one just for fun
    outliers[100, 100] = -hotpixel_intensity

    mock_data = low_order_variability + noise + outliers

    # build an IFUImageModel from these data and give it a wcs
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="LONG")
    hdul[1].data = mock_data

    im = dm.IFUImageModel(hdul)
    ref = create_reference_files(im)
    pipeline = miri.create_pipeline(im, ref)
    wcsobj = wcs.WCS(pipeline)
    im.meta.wcs = wcsobj

    return im


@pytest.fixture(scope="module")
def sci(background):
    """Create science IFUImageMoel for testing.
    Same data as background with different noise.
    """
    hdul = create_hdul(detector="MIRIFULONG", channel="34", band="LONG")
    hdul[1].data = background.data.copy()
    rng = np.random.default_rng(seed=99)
    hdul[1].data += rng.standard_normal(hdul[1].data.shape)
    im = dm.IFUImageModel(hdul)
    im.meta.wcs = background.meta.wcs
    return im


@pytest.fixture(scope="module")
def asn(tmp_cwd_module, background, sci):
    """Create association for testing. Needs at least
    two background images to properly test the step."""

    sci_path = tmp_cwd_module / "sci.fits"
    sci.save(sci_path)

    for root in ["bkg0", "bkg1", "selfcal0", "selfcal1"]:
        path = tmp_cwd_module / (root + ".fits")
        background.save(path)

    asn_table = {
        "asn_pool": "singleton",
        "products": [
            {
                "members": [
                    {"expname": "sci.fits", "exptype": "science"},
                    {"expname": "bkg0.fits", "exptype": "background"},
                    {"expname": "bkg1.fits", "exptype": "background"},
                    {"expname": "selfcal0.fits", "exptype": "selfcal"},
                    {"expname": "selfcal1.fits", "exptype": "selfcal"},
                ]
            }
        ],
    }
    with open(tmp_cwd_module / "tmp_asn.json", "w") as f:
        json.dump(asn_table, f)

    container = dm.open(tmp_cwd_module / "tmp_asn.json")

    return container


def test_input_parsing(asn, sci, background):
    """Test that the input parsing function correctly identifies the science,
    background, and selfcal exposures in the association file
    given association vs imagemodel inputs, and selfcal_list vs not

    Note that science exposure gets added to selfcal_list later, not in parse_inputs
    """

    # basic association case. Both background and selfcal get into the list
    input_sci, selfcal_list, bkg_list = _parse_inputs(asn, [], [])
    assert isinstance(input_sci, dm.IFUImageModel)
    assert len(bkg_list) == 2
    assert len(selfcal_list) == 4

    # association with background_list provided
    input_sci, selfcal_list, bkg_list = _parse_inputs(
        asn,
        [],
        [
            background,
        ]
        * 3,
    )
    assert isinstance(input_sci, dm.IFUImageModel)
    assert len(bkg_list) == 5
    assert len(selfcal_list) == 7

    # association with selfcal_list provided
    input_sci, selfcal_list, bkg_list = _parse_inputs(
        asn,
        [
            background,
        ]
        * 3,
        [],
    )
    assert isinstance(input_sci, dm.IFUImageModel)
    assert len(bkg_list) == 2
    assert len(selfcal_list) == 7

    # single science exposure
    input_sci, selfcal_list, bkg_list = _parse_inputs(sci, [], [])
    assert isinstance(input_sci, dm.IFUImageModel)
    assert len(bkg_list) == 0
    assert len(selfcal_list) == 0

    # single science exposure with selfcal_list and bkg_list provided
    input_sci, selfcal_list, bkg_list = _parse_inputs(
        sci,
        [
            background,
        ]
        * 3,
        [
            background,
        ]
        * 1,
    )
    assert isinstance(input_sci, dm.IFUImageModel)
    assert len(bkg_list) == 1
    assert len(selfcal_list) == 4


def test_background_flagger_mrs(background):
    """
    Ensure the right number of outliers are found, and that
    true outliers are among those.
    """

    bg = background.data

    # pass into the MRSBackgroundFlagger and check it found the right pixels
    flagfrac = 0.001
    result = badpix_selfcal(bg, flagfrac_lower=flagfrac, flagfrac_upper=flagfrac)
    result_tuples = [(i, j) for i, j in zip(*result)]

    # check that the hot pixels were among those flagged
    for idx in outlier_indices:
        assert idx in result_tuples

    # check that the number of flagged pixels is as expected
    assert np.isclose(len(result_tuples) / bg.size, flagfrac * 2, atol=0.0001)


def test_apply_flags(background):
    """
    Ensure that flagged pixels are set to NaN in the data and err arrays,
    and that the DQ flag is set to 1.
    """

    flagged_indices = np.array(outlier_indices).T
    flagged_indices = (flagged_indices[0], flagged_indices[1])

    flagged = apply_flags(background, flagged_indices)

    # check that flagged pixels are NaN in data and err arrays
    for idx in outlier_indices:
        assert np.isnan(flagged.data[idx])
        assert np.isnan(flagged.err[idx])
        assert np.isnan(flagged.var_poisson[idx])
        assert np.isnan(flagged.var_rnoise[idx])
        assert np.isnan(flagged.var_flat[idx])

    # check that DQ flag is set properly
    for idx in outlier_indices:
        assert flagged.dq[idx] == pixel["DO_NOT_USE"] + pixel["OTHER_BAD_PIXEL"]


@pytest.mark.parametrize("dset", ["sci", "asn"])
def test_badpix_selfcal_step(request, dset):
    """Test the badpix_selfcal step. This is a functional test that checks
    that the step runs without error. The output will be checked by the
    regtest.
    """
    input_data = request.getfixturevalue(dset)
    result = BadpixSelfcalStep.call(input_data, skip=False, force_single=True)

    assert result[0].meta.cal_step.badpix_selfcal == "COMPLETE"
    if dset == "sci":
        assert len(result) == 2
        assert isinstance(result[0], dm.IFUImageModel)
        assert len(result[1]) == 0
    else:
        # should return sci, (bkg0, bkg1) but not selfcal0, selfcal1
        assert len(result) == 2
        assert isinstance(result[0], dm.IFUImageModel)
        assert len(result[1]) == 2


def test_expected_fail_sci(sci):
    """Test that the step fails as expected when given only a science exposure
    and force_single is False.
    """
    result = BadpixSelfcalStep.call(sci, skip=False, force_single=False)
    assert result[0].meta.cal_step.badpix_selfcal == "SKIPPED"


@pytest.fixture(scope="module")
def vertical_striping():
    """2-D array with vertical stripes and some outliers, used to test that
    the step is flagging in the correct (along-dispersion) direction
    """
    shp = (102, 96)
    stripes = np.zeros(shp)
    stripes[:, ::4] = hotpixel_intensity

    # add some outliers
    outliers = np.zeros(shp)
    for idx in outlier_indices:
        i, j = int(idx[0] / 10), int(idx[1] / 10)
        outliers[i, j] = hotpixel_intensity
    data = stripes + outliers

    return data


@pytest.mark.parametrize("dispaxis", [None, 1, 2])
def test_dispersion_direction(vertical_striping, dispaxis):
    """
    Test that the flagging operates as expected as a function of direction.

    Of the four outliers, only one lies atop a stripe, so the step should only
    flag one in four as an outlier unless median filter occurs in the along-stripe
    direction.
    """

    flagged_indices = badpix_selfcal(vertical_striping, dispaxis=dispaxis)

    if dispaxis == 2:
        assert len(flagged_indices[0]) == 4

    elif dispaxis == 1:
        assert len(flagged_indices[0]) == 1

    elif dispaxis is None:
        assert len(flagged_indices[0]) == 1
