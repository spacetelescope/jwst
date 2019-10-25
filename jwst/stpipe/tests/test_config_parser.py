import pytest

from .. import config_parser
from ...extern.configobj.configobj import ConfigObj

def test_load_config_file_s3():
    result = config_parser.load_config_file("s3://test-s3-data/pars.asdf")
    assert isinstance(result, ConfigObj)

    with pytest.raises(ValueError):
        config_parser.load_config_file("s3://test-s3-data/missing.asdf")
