import shutil

import pytest

SCRIPTS = [
    "adjust_wcs",
    "asn_edit",
    "asn_gather",
    "asn_make_pool",
    "collect_pipeline_cfgs",
    "create_data",
    "pointing_summary",
    "schemadoc",
    "set_telescope_pointing",
    "set_telescope_pointing.py",
    "set_velocity_aberration",
    "set_velocity_aberration.py",
    "v1_calculate",
    "world_coords",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_installed(script):
    assert shutil.which(script) is not None, f"`{script}` not installed"
