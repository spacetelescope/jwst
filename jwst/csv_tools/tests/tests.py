import os

import pytest

from ..csv_to_hdulist import csv_to_hdulist
from ..csv_to_table import csv_to_table
from ..table_to_json import table_to_json
from . import data


class TestCSVConvert():

    json_mixed_result = r'{"header": {"key": "value", "key2": "value2", "key4": "this is a long one", "key3": "value3"}, "columns": {"first": ["a1", "b1"], "second": ["a2", "b2"], "third": ["a3", "b3"]}}'

    data_path = os.path.split(os.path.abspath(data.__file__))[0]

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_csv_to_hdulist(self):
        """Test basic convserion"""
        path = os.path.join(self.data_path, 'test_csv_mixed_pipe.txt')
        c2h = csv_to_hdulist(path, delimiter='|')
        assert len(c2h) == 2
        assert c2h[0].header['key'] == 'value'
        assert c2h[0].header['key2'] == 'value2'
        assert len(c2h[1].columns) == 3
        assert c2h[1].data['first'][0] == 'a1'

    def test_comment(self):
        """Test with comment"""
        path = os.path.join(self.data_path, 'test_csv_comment.txt')
        c2h = csv_to_hdulist(path, delimiter=',')
        assert c2h[0].header['comkey'] == 'this should have a'
        assert c2h[0].header.comments['comkey'] == 'comment for you'

    @pytest.mark.skip(reason="py3 dicts return results in unpredictable order")
    def test_json(self):
        """Test json"""
        path = os.path.join(self.data_path, 'test_csv_mixed.txt')
        c2json = table_to_json(csv_to_table(path, delimiter=','))
        print(c2json)
        assert c2json == self.json_mixed_result
