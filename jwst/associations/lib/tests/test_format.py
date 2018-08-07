"""FormatTemplate tests"""
from functools import partial
import pytest

from ..format_template import FormatTemplate

@pytest.mark.parametrize(
    'formats, template, expected, fields, errors',
    [
        (['s{:05d}', 's{:s}'], 'astring', 'astring', {}, None),
        (['s{:05d}', 's{:s}'], 'astring', 'astring_s00001', {'field': 1}, None),
        (['s{:05d}', 's{:s}'], 'astring', 'astring_smysource', {'field': "mysource"}, None),
        (['s{:05d}'], 'astring', 'astring_error', {'field': '5.5'}, (RuntimeError,)),
    ]
)
def test_multiple_key_formats(formats, template, expected, fields, errors):
    """Test key formatting options"""

    format_product = FormatTemplate(
        key_formats={
            'field': formats
        }
    )
    statement = partial(format_product, template, **fields)

    if errors is None:
        result = statement()
        assert result == expected
    else:
        with pytest.raises(errors):
            result = statement()
