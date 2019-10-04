"""Test versioning consistency"""

from ..asn_from_list import asn_from_list
from ... import __version__


def test_asn_version():
    """Test version in association is package version"""

    asn = asn_from_list(['a', 'b', 'c'], product_name='aproduct')

    assert asn['code_version'] == __version__
