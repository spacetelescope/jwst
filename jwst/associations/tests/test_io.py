import os
from glob import glob

import pytest
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename

from jwst.associations import load_asn
from jwst.associations.main import Main
from jwst.associations.tests.helpers import combine_pools


@pytest.fixture(scope="module")
def pool():
    """Retrieve pool path"""
    pool_path = get_pkg_data_filename(
        "data/pool_018_all_exptypes.csv", package="jwst.associations.tests"
    )  # Basic pool
    pool = combine_pools(pool_path)

    return pool


@pytest.fixture(scope="module", params=["yaml", "json"])
def make_asns(pool, request, tmp_path_factory):
    asn_format = request.param
    path = str(tmp_path_factory.mktemp("data"))
    if asn_format == "yaml":
        with pytest.warns(
            DeprecationWarning, match="Support for associations as YAML files is deprecated"
        ):
            generated = Main.cli(
                ["-p", path, "-i", "o001", "--save-orphans", "--format", asn_format], pool=pool
            )
    else:
        generated = Main.cli(
            ["-p", path, "-i", "o001", "--save-orphans", "--format", asn_format], pool=pool
        )
    yield generated, path, asn_format


def test_roundtrip(make_asns):
    generated, path, asn_format = make_asns
    asn_files = glob(os.path.join(path, "*." + asn_format))
    assert len(asn_files) == len(generated.associations)

    for asn_file in asn_files:
        with open(asn_file, "r") as asn_fp:
            if asn_format == "yaml":
                with pytest.warns(
                    DeprecationWarning, match="Support for associations as YAML files is deprecated"
                ):
                    load_asn(asn_fp)
            else:
                load_asn(asn_fp)

    orphaned_files = glob(os.path.join(path, "*.csv"))
    assert len(orphaned_files) == 1
    orphaned = Table.read(orphaned_files[0], format="ascii", delimiter="|")
    assert len(orphaned) == len(generated.orphaned)


def test_load_asn_all(make_asns):
    generated, path, asn_format = make_asns
    asn_files = glob(os.path.join(path, "*." + asn_format))
    assert len(asn_files) == len(generated.associations)

    for asn_file in asn_files:
        with open(asn_file, "r") as asn_fp:
            if asn_format == "yaml":
                with pytest.warns(
                    DeprecationWarning, match="Support for associations as YAML files is deprecated"
                ):
                    asns = load_asn(asn_fp, first=False)
            else:
                asns = load_asn(asn_fp, first=False)
        assert len(asns) > 1


def test_warn_invalid_json_as_yaml(tmp_path, log_watcher):
    fname = tmp_path / "bad_asn.json"
    with open(fname, "w") as f:
        # Trailing comma makes this invalid JSON, but YAML accepts it
        f.write(
            '{"asn_rule": "Asn_Lv2Image", "asn_pool": "test", "program": "99999",'
            ' "asn_type": "image2", "products": [{"name": "test_rate",'
            ' "members": [{"expname": "test_rate.fits", "exptype": "science"},]}]}'
        )
    msg = "Association file has json suffix but is invalid JSON"
    watcher = log_watcher(
        "jwst.associations.load_asn",
        message=msg,
        level="warning",
    )
    with open(fname) as f:
        with pytest.warns(DeprecationWarning, match=msg):
            load_asn(f)
    watcher.assert_seen()


def test_warn_force_yaml(tmp_path, log_watcher):
    fname = tmp_path / "force_yaml_asn.json"
    with open(fname, "w") as f:
        # Valid JSON, but we will force load as YAML
        f.write(
            '{"asn_rule": "Asn_Lv2Image", "asn_pool": "test", "program": "99999",'
            ' "asn_type": "image2", "products": [{"name": "test_rate",'
            ' "members": [{"expname": "test_rate.fits", "exptype": "science"}]}]}'
        )
    msg = "Association file has json suffix but is invalid JSON or is force-loaded as YAML"
    watcher = log_watcher(
        "jwst.associations.load_asn",
        message=msg,
        level="warning",
    )
    with open(fname) as f:
        with pytest.warns(DeprecationWarning, match=msg):
            load_asn(f, fmt="yaml")
    watcher.assert_seen()
