"""
Test using SDP-generated pools.
"""

import os
import re
from glob import glob

import pytest

from jwst.associations.lib.diff import compare_asn_files
from jwst.associations.main import Main as asn_generate

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]

# Decompose pool name to retrieve proposal and version id.
pool_regex = re.compile(r"(?P<proposal>jw.+?)_(?P<versionid>.+)_pool")


def parfunc(x):
    """
    For readable param ids.
    ``x`` is tuple of ``(pool, args)``.
    """
    return x[0]


def _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args):
    pool, args = pool_args
    proposal, version_id = pool_regex.match(pool).group("proposal", "versionid")

    input_csv = rtdata.get_data(f"associations/sdp/pools/{pool}.csv")
    rtdata.output = os.curdir  # This test is jailed and parametrized
    rtdata.okify_op = "sdp_pool_copy"  # failures as folder content replacements
    # Create the associations
    with resource_tracker.track(log=request):
        asn_generate.cli(args + ["-p", rtdata.output, "--version-id", version_id, input_csv])
    out_paths = sorted(glob("*.json"))

    # Compare to the truth associations.
    truth_pool_path = f"truth/test_associations_sdp_pools/{pool}/"
    rtdata.truth_remote = truth_pool_path
    truth_paths = sorted(
        [rtdata.get_truth(p) for p in rtdata.data_glob(truth_pool_path, glob="*.json")]
    )
    if truth_paths == []:  # truth dir does not exist
        rtdata.truth_remote = f"{rtdata._inputs_root}/{rtdata.env}/{truth_pool_path}"
    compare_asn_files(out_paths, truth_paths)


# NOTE: These are inflight equivalent approximate replacements for
#       test_associations_standards.py test module.
@pytest.mark.parametrize(
    "pool_args",
    [
        pytest.param(
            ("jw01292_20250316t033413_pool", []), id="pool_006_spec_nirspec_FIXED_SLIT_AND_jw00217"
        ),  # Also test NRS_FSS_VALID_OPTICAL_PATHS
        pytest.param(("jw01467_20250316t025827_pool", []), id="pool_004_wfs"),
        pytest.param(
            ("jw01529_20250316t074500_pool", []),
            id="pool_009_spec_miri_lv2bkg_FIXED_SLIT_AND_pool_011_spec_miri_lv2bkg_lrs",
        ),
        pytest.param(("jw01554_20250316t071400_pool", []), id="imprint_sci_bg_sep_obs_1"),
        pytest.param(("jw01741_20250318t062014_pool", []), id="imprint_n_dithers_obs_6"),
        pytest.param(
            ("jw01794_20250316t024033_pool", []), id="pool_027_nirspec_ifu_nods_NONEWITH_REPEAT"
        ),
        pytest.param(
            ("jw01863_20250316t002754_pool", []), id="pool_010_spec_nirspec_lv2bkg_IFU_BG_DEDICATED"
        ),
        pytest.param(("jw01893_20250316t095836_pool", []), id="imprint_sci_bg_obs_1_2"),
        pytest.param(("jw01958_20250316t041843_pool", []), id="pool_007_spec_miri_MRS"),
        pytest.param(
            ("jw01964_20250316t064614_pool", []),
            id="pool_010_spec_nirspec_lv2bkg_IFU_BG_IN_TARG_GRP",
        ),
        pytest.param(("jw01970_20250317t180031_pool", []), id="pool_017_spec_nirspec_lv2imprint"),
        pytest.param(("jw02113_20250308t170530_pool", []), id="pool_005_spec_niriss"),
        pytest.param(("jw02304_20250316t004555_pool", []), id="pool_026_mir_image_tso"),
        pytest.param(
            ("jw02344_20250318t210812_pool", []), id="pool_027_nirspec_ifu_nods_2NODS_4DTH"
        ),
        pytest.param(("jw02508_20250308t182800_pool", []), id="pool_021_tso_MIRI_LRS"),
        pytest.param(("jw02770_20250409t122318_pool", []), id="imprint_mos_obs_2"),
        pytest.param(("jw02961_20250308t142131_pool", []), id="pool_007_spec_miri_SLITLESS"),
        pytest.param(("jw03596_20250308t161159_pool", []), id="pool_021_tso_NIRISS_SOSS"),
        pytest.param(("jw03777_20250316t024410_pool", []), id="pool_006_spec_nirspec_IFU"),
        pytest.param(("jw03823_20250316t070626_pool", []), id="pool_029_mir_lrsfs_nonod"),
        pytest.param(("jw03969_20250316t131526_pool", []), id="pool_021_tso_NIRSPEC_BRIGHTOBJ"),
        pytest.param(("jw04090_20250316t054542_pool", []), id="pool_013_coron_nircam"),
        pytest.param(("jw04368_20250318t032038_pool", []), id="pool_019_niriss_wfss"),
        pytest.param(("jw04557_20250318t100949_pool", []), id="pool_006_spec_nirspec_MOS"),
        pytest.param(("jw04611_20250308t142406_pool", []), id="pool_014_ami_niriss"),
        pytest.param(
            ("jw05168_20250316t055106_pool", []),
            id="pool_007_spec_miri_FIXED_SLIT_AND_pool_028_mir_lrsfs_nods",
        ),
        pytest.param(("jw05204_20250308t202944_pool", []), id="pool_002_image_miri"),
    ],
    ids=parfunc,
)
def test_std(_jail, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


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
    ids=parfunc,
)
@pytest.mark.slow
def test_sslow(_jail, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw01019_20250316t142947_pool", []),  # FGS_FOCUS + FGS_IMAGE
        ("jw01120_20250317t230449_pool", []),  # NIRSPEC IMAGING + NIRSPEC INTERNAL LAMP
        ("jw01122_uniq_lamp_optical_path_pool", []),  # NRS_FSS_VALID_LAMP_OPTICAL_PATHS
        ("jw01250_20250318t113408_pool", []),  # Moving target in both NRC_IMAGE and MIR_IMAGE
        ("jw01290_20230304t140931_withids_pool", ["-i", "o012", "c1018"]),
        ("jw01411_20250317t143327_pool", ["--DMS"]),  # Coron 5-POINT-SMALL-GRID, DMS for BKG
        ("jw01464_20250318t062856_pool", []),  # WFS&C rules
        ("jw01493_20230307t040130_withids_pool", ["-i", "o003", "c1000"]),
        ("jw01499_20250310t191505_pool", []),  # NIS_DARK
        ("jw01503_20250316t013145_pool", []),  # NIS_LAMP
        ("jw01512_20250316t073810_pool", []),  # NIS_EXTCAL
        ("jw01621_20250318t175957_pool", []),  # NRS_AUTOWAVE + NRS_IFU
        ("jw01678_20240721t195707_pool", []),
        ("jw02064_20230302t112350_withids_pool", ["-i", "o061", "c1008", "c1017"]),
        ("jw02162_20241213t063547_pool", []),
        (
            "jw04225_20241213t150701_pool",
            ["-i", "o001", "o002"],
        ),  # This pair of pools test the DMS flag usage to prevent o-type ASNs when a background c-type candidate is attached to the science exposure.
        ("jw04225_20241213t150701DMS_pool", ["--DMS", "-i", "o001", "o002"]),
        ("jw04462_20250318t100414_pool", []),  # NRS_FSS_VALID_LAMP_OPTICAL_PATHS
        (
            "jw05554_20250528t204800_c1012_pool",
            ["--DMS", "-i", "o009", "o010", "c1012"],
        ),  # This pool checks background behavior with paired MIRI MRS/Imaging exposures
        ("jw04470_20250317t231014_pool", []),  # NIS_IMAGE science program
    ],
    ids=parfunc,
)
def test_sdp(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw01192_o008_pool", ["--DMS", "-i", "o008"]),  # This pool checks imprint behavior.
        ("jw01257_20221201t192226_pool", []),
        ("jw01288_c1005_mostilno12_pool", ["-i", "o003", "c1001", "c1005"]),  # JP-3230
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
@pytest.mark.slow
def test_slow(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)
