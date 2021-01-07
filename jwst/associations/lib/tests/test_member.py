"""Tests for Member"""

from jwst.associations.lib.member import Member


def test_member_from_member():
    """Test member creation from a Member

    Note that this really is a smoke-screen test
    that tests both code and very basic functionality.
    """
    data = {'exp_name': 'fred', 'exp_type': 'not_science'}
    item = 'This is an item'
    member = Member(data, item=item)
    dup = Member(member)
    assert dup.data == member.data
    assert dup.item == member.item
