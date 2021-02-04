import pytest

from stpipe import config_parser
from stpipe.extern.configobj.configobj import ConfigObj

from jwst.stpipe.tests.steps import MakeListStep


def test_load_config_file_s3(s3_root_dir):
    path = str(s3_root_dir.join("pars-makeliststep.asdf"))
    MakeListStep(par1=42.0, par2="foo").export_config(path)

    result = config_parser.load_config_file("s3://test-s3-data/pars-makeliststep.asdf")
    assert isinstance(result, ConfigObj)

    with pytest.raises(ValueError):
        config_parser.load_config_file("s3://test-s3-data/missing.asdf")
