from collections import MutableMapping

from ..lib.dictwithattrs import DictWithAttributes
from ..lib.member import Member

import pytest

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
        pytest.raises(KeyError, lambda x: x['a'], da)
        pytest.raises(AttributeError, lambda x: x.a, da)
