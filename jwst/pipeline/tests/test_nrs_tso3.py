"""Test calwebb_tso3 with NRS"""

from collections import defaultdict
from glob import glob
from os import path

from .helpers import (
    abspath,
    mk_tmp_dirs,
    require_bigdata,
    runslow,
    update_asn_basedir,
)

from ...associations import load_asn
from ...stpipe.step import (Step, remove_suffix)

DATAPATH = abspath(
    path.join('$TEST_BIGDATA', 'pipelines', 'nirspec_caltso3')
)


@runslow
@require_bigdata
def test_run_full(mk_tmp_dirs):
    """Test a full run"""
    tmp_current_path, tmp_data_path, tmp_config_path = mk_tmp_dirs

    asn_path = update_asn_basedir(
        path.join(DATAPATH, 'jw93056-o001_20180519t045329_tso3_001_asn.json'),
        root=DATAPATH
    )
    args = [
        path.join(DATAPATH, 'cfgs', 'calwebb_tso3.cfg'),
        asn_path,
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
    product_name_file = product_name + '_whtlt.ecsv'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    product_name_file = product_name + '_x1dints.fits'
    assert product_name_file in output_files
    output_files.remove(product_name_file)

    # Check Level2 products
    for member in members_by_type['science']:
        name, ext = path.splitext(path.split(member)[1])
        name, separator = remove_suffix(name)
        name = name + separator + acid + separator + 'crfints' + ext
        assert name in output_files
        output_files.remove(name)

    # If there are files left, this is an error
    assert len(output_files) == 0
