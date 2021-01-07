import os
import sys

import pytest

from ..csv_to_hdulist import csv_to_hdulist
from ..csv_to_table import csv_to_table
from ..table_to_json import table_to_json
from . import data


json_mixed_result = "".join([
    r'{"header": {"key": "value", "key2": "value2", '
    r'"key4": "this is a long one", "key3": "value3"}, "columns": '
    r'{"first": ["a1", "b1"], "second": ["a2", "b2"], "third": ["a3", "b3"]}}'
    ])

data_path = os.path.split(os.path.abspath(data.__file__))[0]


@pytest.mark.xfail(reason="csv_to_hdulist is broken")
def test_csv_to_hdulist():
    """Test basic convserion"""
    path = os.path.join(data_path, 'test_csv_mixed_pipe.txt')
    c2h = csv_to_hdulist(path, delimiter='|')
    assert len(c2h) == 2
    assert c2h[0].header['key'] == 'value'
    assert c2h[0].header['key2'] == 'value2'
    assert len(c2h[1].columns) == 3
    assert c2h[1].data['first'][0] == 'a1'

@pytest.mark.xfail(reason="csv_to_hdulist is broken")
def test_comment():
    """Test with comment"""
    path = os.path.join(data_path, 'test_csv_comment.txt')
    c2h = csv_to_hdulist(path, delimiter=',')
    assert c2h[0].header['comkey'] == 'this should have a'
    assert c2h[0].header.comments['comkey'] == 'comment for you'

@pytest.mark.skipif(sys.version_info < (3,6),
                    reason="requires python>=3.6 so dicts are ordered")
def test_json():
    """Test json"""
    path = os.path.join(data_path, 'test_csv_mixed.txt')
    c2json = table_to_json(csv_to_table(path, delimiter=','))
    assert c2json == json_mixed_result
