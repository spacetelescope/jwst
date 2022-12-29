"""Test against Level2 standard associations
"""
import os
from pathlib import Path
import pytest

from jwst.associations.lib.diff import (
    compare_asn_files,
)
from jwst.associations.tests.helpers import (
    combine_pools,
    t_path
)
from jwst.lib.file_utils import pushdir

from jwst.associations.main import Main

# #################
# Setup environment
# #################

# Produce Level2b only associations
LV2_ONLY_ARGS = [
    '-r',
    t_path('../lib/rules_level2b.py'),
    '--ignore-default',
]

# Produce Level3 only associations
LV3_ONLY_ARGS = [
    '-r',
    t_path('../lib/rules_level3.py'),
    '--ignore-default',
]

# Produce general associations
DEF_ARGS = []

# Define the standards


class MakePars():
    """Setup the test parameters """

    def __init__(
            self,
            pool_root,
            main_args=DEF_ARGS,
            source=None,
            outdir=None,
            execute=True,
            xfail=None,
            slow=False,
    ):
        self.pool_root = pool_root
        self.main_args = main_args
        self.source = source
        self.outdir = outdir
        self.execute = execute
        self.xfail = xfail
        self.slow = slow


standards = [
    MakePars('pool_002_image_miri', main_args=LV3_ONLY_ARGS),
    MakePars('pool_004_wfs', slow=True),
    MakePars('pool_005_spec_niriss'),
    MakePars('pool_006_spec_nirspec'),
    MakePars('pool_007_spec_miri'),
    MakePars('pool_009_spec_miri_lv2bkg'),
    MakePars('pool_010_spec_nirspec_lv2bkg'),
    MakePars('pool_011_spec_miri_lv2bkg_lrs'),
    MakePars('pool_013_coron_nircam', slow=True),
    MakePars('pool_014_ami_niriss'),
    MakePars('pool_015_spec_nirspec_lv2bkg_reversed', main_args=LV2_ONLY_ARGS),
    MakePars('pool_016_spec_nirspec_lv2bkg_double', main_args=LV2_ONLY_ARGS),
    MakePars('pool_017_spec_nirspec_lv2imprint'),
    MakePars('pool_018_all_exptypes', main_args=LV2_ONLY_ARGS),
    MakePars('pool_019_niriss_wfss'),
    MakePars('pool_020_00009_image_miri'),
    MakePars('pool_021_tso'),
    MakePars('pool_022_tso_noflag'),
    MakePars('pool_023_nirspec_msa_3nod', main_args=LV2_ONLY_ARGS),
    MakePars('pool_024_nirspec_fss_nods'),
    MakePars('pool_025_nirspec_fss_nod_chop'),
    MakePars('pool_026_mir_image_tso'),
    MakePars('pool_027_nirspec_ifu_nods'),
    MakePars('pool_028_mir_lrsfs_nods'),
    MakePars('pool_029_mir_lrsfs_nonod'),
    MakePars('pool_030_mir_lrs_nods_bkg'),
    MakePars('pool_031_mir_lrs_nonod_bkg'),
    MakePars('pool_032_nircam_wfss'),
]


# #####
# Tests
# #####
def generate_id(value):
    """Generate test ids based on the parametrized input"""
    return value.pool_root


input_loc = 'associations'
test_dir = 'standards'
ref_loc = [input_loc, test_dir, 'truth']


@pytest.mark.bigdata
@pytest.mark.filterwarnings('error')
@pytest.mark.parametrize('standard_pars', standards, ids=generate_id)
def test_against_standard(rtdata, standard_pars, slow):
    """Compare a generated association against a standard
    Success is when no other AssertionError occurs.
    """
    if standard_pars.xfail is not None:
        pytest.xfail(reason=standard_pars.xfail)

    if standard_pars.slow and not slow:
        pytest.skip(f'Pool {standard_pars.pool_root} requires "--slow" option')

    # Tell rtdata to handle failures as folder replacements.
    rtdata.okify_op = 'folder_copy'

    # Start test
    pool_root = standard_pars.pool_root
    cwd = Path(pool_root)
    cwd.mkdir()
    with pushdir(cwd):

        # Create the associations
        output_path = Path(pool_root)
        output_path.mkdir()
        rtdata.output = str(output_path)
        version_id = standard_pars.pool_root.replace('_', '-')
        args = standard_pars.main_args + [
            '-p', str(output_path),
            '--version-id', version_id,
        ]
        pool = combine_pools([
            t_path(Path('data') / (standard_pars.pool_root + '.csv'))
        ])
        Main.cli(args, pool=pool)

        # Retrieve the truth files
        truth_paths = []
        truth_pool_path = '/'.join(ref_loc) + '/' + standard_pars.pool_root
        for path in rtdata.data_glob(truth_pool_path, glob='*.json'):
            truth_paths.append(rtdata.get_truth(path))
        rtdata.truth_remote = os.path.join(rtdata._inputs_root, rtdata.env, truth_pool_path)

        # Compare the association sets.
        try:
            compare_asn_files(output_path.glob('*.json'), truth_paths)
        except AssertionError:
            if standard_pars.xfail:
                pytest.xfail(standard_pars.xfail)
            else:
                raise
