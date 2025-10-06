"""Test using SDP-generated pools."""

import pytest

from jwst.regtest.associations_sdp_pools.assoc_rt_helpers import assoc_sdp_against_standard, parfunc

# Mark all tests in this module
pytestmark = [pytest.mark.bigdata]


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
)
def test_std(_jail, rtdata, resource_tracker, request, pool_args):
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)


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
        ("jw04470_20250317t231014_pool", []),  # NIS_IMAGE science program
        (
            "jw05554_20250528t204800_c1012_pool",
            ["--DMS", "-i", "o009", "o010", "c1012"],
        ),  # This pool checks background behavior with paired MIRI MRS/Imaging exposures
    ],
    ids=parfunc,
)
def test_sdp(tmp_cwd, rtdata, resource_tracker, request, pool_args):
    assoc_sdp_against_standard(rtdata, resource_tracker, request, pool_args)
