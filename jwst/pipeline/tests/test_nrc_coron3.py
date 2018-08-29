"""Test calwebb_coron3 with NIRCam"""

from collections import defaultdict
from glob import glob
from os import path

from .helpers import (
    SCRIPT_PATH,
    SCRIPT_DATA_PATH,
    abspath,
    mk_tmp_dirs,
    require_bigdata,
    runslow,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import (Step, remove_suffix)

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'nircam', 'test_coron3')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw99999-a3001_20170327t121212_coron3_001_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(SCRIPT_DATA_PATH, 'cfgs', 'calwebb_coron3.cfg'),
        asn_path,
        '--steps.align_refs.override_psfmask=' + path.join(DATAPATH, 'jwst_nircam_psfmask_somb.fits'),
        '--steps.resample.skip=true'
    ]

    Step.from_cmdline(args)

    # Now test for file existence. Get the association
    with open(asn_path) as fh:
        asn = load_asn(fh)
    acid = asn['asn_id']
    product = asn['products'][0]
    product_name = product['name']
    members_by_type = defaultdict(list)
    for member in product['members']:
        expname = path.split(member['expname'])[1]
        members_by_type[member['exptype'].lower()].append(expname)

    output_files = glob('*')
    print('Created files ares: {}'.format(output_files))

    # Check Level3 products
    product_name_file = product_name + '_psfstack.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    product_name_file = product_name + '_i2d.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    # Check Level2 products
    for member in members_by_type['science']:
        basename, ext = path.splitext(path.split(member)[1])
        basename, separator = remove_suffix(basename)

        name = basename + separator + acid + separator + 'psfalign' + ext
        assert name in output_files
        output_files.remove(name)

        name = basename + separator + acid + separator + 'psfsub' + ext
        assert name in output_files
        output_files.remove(name)

        name = basename + separator + acid + separator + 'crfints' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0
