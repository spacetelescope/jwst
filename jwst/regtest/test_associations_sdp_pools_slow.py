"""Test using SDP-generated pools (slow)."""

import pytest

from jwst.regtest.assoc_rt_helpers import assoc_sdp_against_standard, parfunc

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.slow]


# NOTE: These are inflight equivalent approximate replacements for
#       test_associations_standards.py test module (slow).
@pytest.mark.parametrize(
    "pool_args",
    [
        pytest.param(("jw01523_20250321t155408_pool", []), id="pool_009_spec_miri_lv2bkg_MRS"),
        pytest.param(("jw01928_20250318t105523_pool", []), id="imprint_single_obs_1"),
        pytest.param(
            ("jw02084_20250320t081615_pool", []), id="pool_021_tso_NIRCAM_TSIMAGE_TSGRISM"
        ),
        pytest.param(("jw03383_20250307t235057_pool", []), id="pool_034_wfss_parallel_NIRISS"),
        pytest.param(
            ("jw03522_20250318t192624_pool", []), id="pool_024_nirspec_fss_nods_AND_SUBPXPTS4"
        ),
        pytest.param(("jw03702_20250320t002824_pool", []), id="imprint_2n_dithers_obs_1_2"),
        pytest.param(("jw05221_20250411t182730_pool", []), id="pool_032_nircam_wfss"),
        pytest.param(("jw05308_20250319t215735_pool", []), id="pool_027_nirspec_ifu_nods_4NODS"),
        pytest.param(
            ("jw05398_20250312t000159_default_pool", ["-i", "o041", "o042", "o043", "c1016"]),
            id="pool_034_wfss_parallel_NIRCAM_DEFAULT",
        ),  # This is "default" scenario: row grism, column grism, and direct image.
        pytest.param(
            (
                "jw05398_20250312t000159_3r2d_pool",
                ["-i", "o030", "o031", "o032", "o033", "o035", "c1012", "c1013", "c1014"],
            ),
            id="pool_034_wfss_parallel_NIRCAM_3ROW_2DIRECT",
        ),  # Three row grism images and two direct images.
        pytest.param(
            ("jw05398_20250312t000159_row_pool", ["-i", "o036", "o037", "c1029"]),
            id="pool_034_wfss_parallel_NIRCAM_ROW_ONLY",
        ),  # Only row grism and direct image, but no column grism image.
    ],
)
def test_sslow(_jail, rtdata, resource_tracker, request, pool_args):
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw01192_o008_pool", ["--DMS", "-i", "o008"]),  # This pool checks imprint behavior.
        ("jw01257_20221201t192226_pool", []),
        ("jw01288_2025_c1005_mostilno12_pool", ["-i", "o003", "c1001", "c1005"]),  # JP-3230
        ("jw01290_20230304t140931_pool", []),
        ("jw01480_20250319t173819_pool", []),  # NRC_GRISM, NRC_TSIMAGE, NRC_TSGRISM
        ("jw01493_20230307t040130_pool", []),
        ("jw02064_20230302t112350_pool", []),
        (
            "jw02739_20230710t150016_pool",
            ["-i", "c1000"],
        ),  # This association tests the Asn_Lv3ImageMosaic rule
        ("jw02739_20230710t150016_o010_pool", ["-i", "o010"]),  # NRC_IMAGE
        ("jw03855_20241103t042455_pool", []),
        ("jw04237_20250321t192812_pool", []),  # Lvl 3 MIRI MRS BKG
        ("jw04237_20250321t192812_dms_pool", ["--DMS"]),  # Lvl 3 MIRI MRS BKG with DMS
    ],
    ids=parfunc,
)
def test_slow(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)
