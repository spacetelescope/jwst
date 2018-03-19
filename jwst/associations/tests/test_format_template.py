"""Tests of the format_template module"""

from ..lib.format_template import FormatTemplate


def test_basics():

    # Just defaults
    template = 'name="{name}" value="{value}"'
    fmt = FormatTemplate()
    result = fmt(template)
    assert result == 'name="{name}" value="{value}"'

    # With values given
    result = fmt(template, name='fred', value='great')
    assert result == 'name="fred" value="great"'

    # But wait, too many values given:
    result = fmt(template, name='fred', value='great', extra='more')
    assert result == 'name="fred" value="great"_more'

    # And with too many and not enough:
    result = fmt(template, value='great', extra='more')
    assert result == 'name="{name}" value="great"_more'

    # Nothing should be added if value is None
    result = fmt(template, value=None)
    assert result == 'name="{name}" value=""'

    # With a different separator:
    fmt.separator = '---'
    result = fmt(template, name='fred', value='great', extra='more')
    assert result == 'name="fred" value="great"---more'

    # Initializing with a different separator:
    fmt_newsep = FormatTemplate(separator='_now-with_')
    result = fmt_newsep(template, name='fred', value='great', extra='more')
    assert result == 'name="fred" value="great"_now-with_more'

    # Setup preformatting
    key_formats = {'value': 'pre_{:s}_format'}
    fmt_preformat = FormatTemplate(key_formats=key_formats)
    result = fmt_preformat(template, name='fred', value='great')
    assert result == 'name="fred" value="pre_great_format"'

    # Ignore unused replacements
    fmt = FormatTemplate(remove_unused=True)
    result = fmt(template)
    assert result == 'name="" value=""'
