import shutil

import pytest

SCRIPTS = [
    'adjust_wcs',
    'asn_edit',
    'asn_gather',
    'asn_make_pool',
    'assign_wcs',
    'collect_pipeline_cfgs',
    'coron',
    'create_data',
    'cube_build',
    'dark_current',
    'data_generate',
    'dqinit',
    'flatfieldcorr',
    'fringecorr',
    'ipc',
    'jump',
    'linearitycorr',
    'make_header',
    'migrate_data',
    'minimum_deps',
    'move_wcs',
    'okify_regtests',
    'outlier_detection',
    'persistencecorr',
    'photomcorr',
    'pointing_summary',
    'rampfitcorr',
    'refpix',
    'resample',
    'saturationcorr',
    'schema_editor',
    'schemadoc',
    'set_telescope_pointing',
    'set_telescope_pointing.py',
    'set_velocity_aberration',
    'set_velocity_aberration.py',
    'straylight',
    'superbias',
    'v1_calculate',
    'verify_install_requires',
    'world_coords',
]


@pytest.mark.parametrize('script', SCRIPTS)
def test_script_installed(script):
    assert shutil.which(script) is not None, f'`{script}` not installed'
