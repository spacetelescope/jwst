"""
Test using SDP-generated pools.
"""

import os
import re
from glob import glob

import pytest

from jwst.associations.lib.diff import MultiDiffError, compare_asn_files
from jwst.associations.main import Main as asn_generate

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata, pytest.mark.filterwarnings("error")]

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
        pytest.param(("jw01523_20250321t155408_pool", []), id="pool_009_spec_miri_lv2bkg_MRS"),
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
        pytest.param(("jw01928_20250318t105523_pool", []), id="imprint_single_obs_1"),
        pytest.param(("jw01958_20250316t041843_pool", []), id="pool_007_spec_miri_MRS"),
        pytest.param(
            ("jw01964_20250316t064614_pool", []),
            id="pool_010_spec_nirspec_lv2bkg_IFU_BG_IN_TARG_GRP",
        ),
        pytest.param(("jw01970_20250317t180031_pool", []), id="pool_017_spec_nirspec_lv2imprint"),
        pytest.param(
            ("jw02084_20250320t081615_pool", []), id="pool_021_tso_NIRCAM_TSIMAGE_TSGRISM"
        ),
        pytest.param(("jw02113_20250308t170530_pool", []), id="pool_005_spec_niriss"),
        pytest.param(("jw02304_20250316t004555_pool", []), id="pool_026_mir_image_tso"),
        pytest.param(
            ("jw02344_20250318t210812_pool", []), id="pool_027_nirspec_ifu_nods_2NODS_4DTH"
        ),
        pytest.param(("jw02508_20250308t182800_pool", []), id="pool_021_tso_MIRI_LRS"),
        pytest.param(("jw02770_20250409t122318_pool", []), id="imprint_mos_obs_2"),
        pytest.param(("jw02961_20250308t142131_pool", []), id="pool_007_spec_miri_SLITLESS"),
        pytest.param(("jw03522_20250318t192624_pool", []), id="pool_024_nirspec_fss_nods"),
        pytest.param(("jw03596_20250308t161159_pool", []), id="pool_021_tso_NIRISS_SOSS"),
        pytest.param(("jw03702_20250320t002824_pool", []), id="imprint_2n_dithers_obs_1_2"),
        pytest.param(("jw03777_20250316t024410_pool", []), id="pool_006_spec_nirspec_IFU"),
        pytest.param(("jw03823_20250316t070626_pool", []), id="pool_029_mir_lrsfs_nonod"),
        pytest.param(("jw03969_20250316t131526_pool", []), id="pool_021_tso_NIRSPEC_BRIGHTOBJ"),
        pytest.param(("jw04368_20250318t032038_pool", []), id="pool_019_niriss_wfss"),
        pytest.param(("jw04557_20250318t100949_pool", []), id="pool_006_spec_nirspec_MOS"),
        pytest.param(("jw04611_20250308t142406_pool", []), id="pool_014_ami_niriss"),
        pytest.param(
            ("jw05168_20250316t055106_pool", []),
            id="pool_007_spec_miri_FIXED_SLIT_AND_pool_028_mir_lrsfs_nods",
        ),
        pytest.param(("jw05204_20250308t202944_pool", []), id="pool_002_image_miri"),
        pytest.param(("jw05221_20250411t182730_pool", []), id="pool_032_nircam_wfss"),
        pytest.param(("jw05308_20250319t215735_pool", []), id="pool_027_nirspec_ifu_nods_4NODS"),
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
        pytest.param(("jw01467_20250316t025827_pool", []), id="pool_004_wfs"),
        pytest.param(("jw03383_20250307t235057_pool", []), id="pool_034_wfss_parallel_NIRISS"),
        pytest.param(("jw04090_20250316t054542_pool", []), id="pool_013_coron_nircam"),
        pytest.param(
            ("jw05398_20250312t000159_pool", ["-i", "o041", "o042", "o043", "c1016"]),
            id="pool_034_wfss_parallel_NIRCAM_DEFAULT",
        ),  # This is "default" scenario: row grism, column grism, and direct image.
        pytest.param(
            ("jw05398_20250312t000159_pool", ["-i", "o036", "o037", "c1029"]),
            id="pool_034_wfss_parallel_NIRCAM_ROW_ONLY",
        ),  # Only row grism and direct image, but no column grism image.
        pytest.param(
            (
                "jw05398_20250312t000159_pool",
                ["-i", "o030", "o031", "o032", "o033", "o035", "c1012", "c1013", "c1014"],
            ),
            id="pool_034_wfss_parallel_NIRCAM_3ROW_2DIRECT",
        ),  # Three row grism images and two direct images.
    ],
    ids=parfunc,
)
@pytest.mark.slow
def test_sslow(_jail, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw00623_20200918t091537_o055_pool", []),
        ("jw00624_20190205t031003_pool", []),
        ("jw00625_20190603t233254_pool", []),
        ("jw00626_20190605t025021_pool", []),
        ("jw00632_20210921t193255_pool", []),
        ("jw00660_20210530t152401_pool", []),
        ("jw00663_20221218t111937_pool", ["-i", "o004", "c1000"]),
        ("jw00676_20210403t114320_c1007_pool", ["--DMS"]),
        ("jw00711_20181027T043250_pool", []),
        ("jw00791_20181019T221608_pool", []),
        ("jw00809_20220124t091748_pool", []),
        ("jw00818_20230407t030411_pool", []),
        ("jw00839_20221220t025418_pool", ["-i", "o002", "c1000"]),
        ("jw01122_uniq_lamp_optical_path_pool", []),  # NRS_FSS_VALID_LAMP_OPTICAL_PATHS
        ("jw01192_o008_pool", ["--DMS", "-i", "o008"]),  # This pool checks imprint behavior.
        ("jw01290_20230304t140931_withids_pool", ["-i", "o012", "c1018"]),
        ("jw01493_20230307t040130_pool", []),
        ("jw01493_20230307t040130_withids_pool", ["-i", "o003", "c1000"]),
        ("jw01678_20240721t195707_pool", []),
        ("jw02064_20230302t112350_withids_pool", ["-i", "o061", "c1008", "c1017"]),
        ("jw02162_20241213t063547_pool", []),
        (
            "jw02739_20230710t150016_pool",
            ["-i", "c1000"],
        ),  # This association tests the Asn_Lv3ImageMosaic rule
        ("jw03855_20241103t042455_pool", []),
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
        ("jw10002_20171108T041457_pool", []),
        ("jw10004_20171108T090028_pool", []),
        ("jw10005_20181020T033546_pool", []),
        ("jw84600_20180824T212338_pool", []),
        ("jw84600_20180824T212338-valid-msametfl_pool", []),
        ("jw84700_subpxpts_pool", []),
        ("jw86600_20171108T043532_pool", []),
        ("jw87600_20180824T213416_pool", []),
        ("jw87800_20180824T214549_pool", []),
        ("jw88600_20171108T051731_pool", []),
        ("jw90001_20171108T043223_pool", []),
        ("jw90002_20171108T051832_pool", []),
        ("jw90003_20171108T043229_pool", []),
        ("jw93025_20171108T062313_pool", []),
        ("jw93045_20171108T045547_pool", []),
        ("jw93056_20171108T060152_pool", []),
        ("jw93125_20171108T045251_pool", []),
        ("jw93125_20171108T045251-valid-msametfl_pool", []),
        ("jw95065_20171108T043122_pool", []),
        ("jw95065_20171108T043122-valid-msametfl_pool", []),
        ("jw95115_20171108T041653_pool", []),
        ("jw95175_20171108T044201_pool", []),
        ("jw96090_20171108T041421_pool", []),
        ("jw96213_20171108T053001_pool", []),
        ("jw96215_20180602T170215_pool", []),
        ("jw96691_20171108T060643_pool", []),
        ("jw97012_20171108T064939_pool", []),
        ("jw97013_20171108T044023_pool", []),
        ("jw98005_20171108T041409_pool", []),
    ],
    ids=parfunc,
)
def test_sdp(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw00623_20190607t021101_pool", []),
        ("jw00628_20191102t153956_pool", []),
        ("jw00629_20190605t025157_pool", []),
        ("jw00676_20210403t114320_pool", []),
        ("jw01194_20230115t113819_pool", ["--DMS"]),  # This pool checks coronagraphy associations
        ("jw01257_20221201t192226_pool", []),
        ("jw01290_20230304t140931_pool", []),
        ("jw02064_20230302t112350_pool", []),
        ("jw82600_20180921T023255_pool", []),
        ("jw93065_20171108T041402_pool", []),
        ("jw93135_20171108T041617_pool", []),
        ("jw93135_20171108T041617-fixed_pool", []),
        ("jw94015_20171108T041516_pool", []),
    ],
    ids=parfunc,
)
@pytest.mark.slow
def test_slow(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw80600_20171108T041522_pool", []),  # PR 3450
        ("jw98010_20171108T062332_pool", []),  # PR 3450
    ],
    ids=parfunc,
)
def test_fail(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    with pytest.raises(MultiDiffError):
        _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


@pytest.mark.parametrize(
    "pool_args",
    [
        ("jw01288_c1005_mostilno12_pool", ["-i", "o003", "c1001", "c1005"]),  # JP-3230
    ],
    ids=parfunc,
)
@pytest.mark.slow
def test_fslow(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    with pytest.raises(MultiDiffError):
        _assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)
