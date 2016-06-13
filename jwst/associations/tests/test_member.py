from nose.tools import assert_raises

from collections import MutableMapping

from jwst.associations.lib.dictwithattrs import DictWithAttributes
from jwst.associations.lib.member import Member

class TestMember():
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_dictwithattributes(self):
        da = DictWithAttributes()
        assert len(da) == 0
        da.a = 1
        assert da.a == 1
        assert da['a'] == 1
        da['b'] = 2
        assert da.b == 2
        assert da['b'] == 2
        assert len(da) == 2
        assert isinstance(da, MutableMapping)

    def test_member_nomap(self):
        da = Member()
        assert len(da) == 0
        da.a = 1
        assert da.a == 1
        assert da['a'] == 1
        da['b'] = 2
        assert da.b == 2
        assert da['b'] == 2
        assert len(da) == 2
        assert isinstance(da, MutableMapping)

    def test_member_map(self):
        indict = {'a': 1, 'b': 2, 'c': 3}
        map = {'a': 'aa', 'b': 'bb'}
        da = Member(indict, map)
        assert da.aa == 1
        assert da.c == 3
        assert da['aa'] == 1
        assert da['c'] == 3
        assert_raises(KeyError, lambda x: x['a'], da)
        assert_raises(AttributeError, lambda x: x.a, da)
