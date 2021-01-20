"""FormatTemplate tests"""
from functools import partial
import pytest

from ..format_template import FormatTemplate

@pytest.mark.parametrize(
    'key_formats, template, expected, fields, errors',
    [
        # No replacement at all.
        (
            None,
            'name="{name}" value="{value}"',
            'name="{name}" value="{value}"',
            {},
            None
        ),

        # Basic replacement
        (
            None,
            'name="{name}" value="{value}"',
            'name="fred" value="great"',
            {'name': 'fred', 'value': 'great'},
            None
        ),

        # But wait, too many values given:
        (
            None,
            'name="{name}" value="{value}"',
            'name="fred" value="great"_more',
            {'name': 'fred', 'value': 'great', 'extra': 'more'},
            None
        ),

        # And with too many and not enough:
        (
            None,
            'name="{name}" value="{value}"',
            'name="{name}" value="great"_more',
            {'value': 'great', 'extra': 'more'},
            None
        ),

        # Nothing should be added if value is None
        (
            None,
            'name="{name}" value="{value}"',
            'name="{name}" value=""',
            {'value': None},
            None
        ),

        # Multiple key formats, using none.
        (
            {
                'field': ['s{:05d}', 's{:s}'],
            },
            'astring',
            'astring',
            {},
            None
        ),

        # Multiple key formats, using first.
        (
            {
                'field': ['s{:05d}', 's{:s}'],
            },
            'astring',
            'astring_s00001',
            {'field': 1},
            None
        ),

        # Multiple key formats, using second.
        (
            {
                'field': ['s{:05d}', 's{:s}']
            },
            'astring',
            'astring_smysource',
            {'field': "mysource"},
            None
        ),

        # No matching formats is an error.
        (
            {
                'field': ['s{:05d}']
            },
            'astring',
            'astring_error',
            {'field': '5.5'},
            (RuntimeError,)
        ),
    ]
)
def test_basics(key_formats, template, expected, fields, errors):
    """Test all basic formatting options"""

    format_product = FormatTemplate(key_formats=key_formats)
    statement = partial(format_product, template, **fields)

    if errors is None:
        result = statement()
        assert result == expected
    else:
        with pytest.raises(errors):
            result = statement()


def test_separators():
    """Test changing separators"""
    template = 'name="{name}" value="{value}"'
    fmt = FormatTemplate()

    # With a different separator:
    fmt.separator = '---'
    result = fmt(template, name='fred', value='great', extra='more')
    assert result == 'name="fred" value="great"---more'

    # Initializing with a different separator:
    fmt_newsep = FormatTemplate(separator='_now-with_')
    result = fmt_newsep(template, name='fred', value='great', extra='more')
    assert result == 'name="fred" value="great"_now-with_more'


def test_allow_unknown():
    """Keep None values in string"""
    template = 'name="{name}" value="{value}"'
    fmt = FormatTemplate(remove_unused=False)
    result = fmt(template)
    assert result == template
