import nose.tools as nt
from nose import SkipTest

from jwst.csv_tools.csv_to_hdulist import *
from jwst.csv_tools.table_to_json import table_to_json

class TestCSVConvert():

    json_mixed_result = r'{"header": {"key": "value", "key2": "value2", "key4": "this is a long one", "key3": "value3"}, "columns": {"first": ["a1", "b1"], "second": ["a2", "b2"], "third": ["a3", "b3"]}}'

    def setUp(self):
        pass

    def tearDown(self):
        pass

    # Basic convserion
    def test_csv_to_hdulist(self):
        c2h = csv_to_hdulist('tests/data/test_csv_mixed_pipe.txt', delimiter='|')
        assert len(c2h) == 2
        assert c2h[0].header['key'] == 'value'
        assert c2h[0].header['key2'] == 'value2'
        assert len(c2h[1].columns) == 3
        assert c2h[1].data['first'][0] == 'a1'

    # Test with comment
    def test_comment(self):
        c2h = csv_to_hdulist('tests/data/test_csv_comment.txt', delimiter=',')
        assert c2h[0].header['comkey'] == 'this should have a'
        assert c2h[0].header.comments['comkey'] == 'comment for you'

    # Test json
    def test_json(self):
        c2json = table_to_json(csv_to_table('tests/data/test_csv_mixed.txt', delimiter=','))
        assert c2json == self.json_mixed_result

# Utility tests.
def check_in_list(element, alist):
    assert element in alist

def check_equal(left, right):
    assert left == right
    
