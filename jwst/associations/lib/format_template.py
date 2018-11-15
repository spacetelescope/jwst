"""FormatTemplate

Format template string allowing partial formatiing.
"""
from collections import defaultdict
from string import Formatter

__all__ = ['FormatTemplate']


# Define conversion based on format type
CONVERSION = {
    '%': float,
    'b': int,
    'c': int,
    'd': int,
    'e': float,
    'E': float,
    'f': float,
    'F': float,
    'g': float,
    'G': float,
    'n': int,
    'o': int,
    's': str,
    'x': int,
    'X': int,
}


class FormatTemplate(Formatter):
    """Format a template

    Parameters
    ----------
    template : str
        A Python format string

    separator : string
        Separater to use for values which have no
        matching replacement strings

    key_formats : dict or None
        A dict of key-specific formatting where the value will
        be pre-formatted before being passed to the final format
        string. Each format will be tried until success.

    remove_unused : bool
        By default, unused replacement fields are left in the
        result, for use in subsequent replacement usage.
        If True, such fields are removed from the result.

    kwargs : dict or named parameters
        The key/values pairs to fill into the Python format string

    Returns
    -------
    str
        The formatted string

    Notes
    -----
    This differences from Pythons `format` method are:
        - If a replacement field does not have a given value,
          the replacement field is left in the result
        - If a key/value pair is present but has no replacement field,
          the value is simply appended.
        - Template can only use named replacement fields.

    Examples
    --------
    The basic example:
    >>> template = 'name="{name}" value="{value}"'
    >>> fmt = FormatTemplate()
    >>> fmt(template)
    'name="{name}" value="{value}"'

    But with actual values given:
    >>> fmt(template, name='fred', value='great')
    'name="fred" value="great"'

    But wait, too many values given:
    >>> fmt(template, name='fred', value='great', extra='more')
    'name="fred" value="great"_more'

    And with too many and not enough:
    >>> fmt(template, value='great', extra='more')
    'name="{name}" value="great"_more'

    With a different separator:
    >>> fmt.separator = '---'
    >>> fmt(template, name='fred', value='great', extra='more')
    'name="fred" value="great"---more'

    Initializing with a different separator:
    >>> fmt_newsep = FormatTemplate(separator='_now-with_')
    >>> fmt_newsep(template, name='fred', value='great', extra='more')
    'name="fred" value="great"_now-with_more'

    Setup preformatting
    >>> key_formats = {'value': ['pre_{:s}_format']}
    >>> fmt_preformat = FormatTemplate(key_formats=key_formats)
    >>> fmt_preformat(template, name='fred', value='great')
    'name="fred" value="pre_great_format"'
    """
    def __init__(self, separator='_', key_formats=None, remove_unused=False):
        """Inialize class

        Parameters
        ----------
        separator : str
            For key/value pairs given that do not have a
            replacement field, the values are appened to
            the string using this separator.

        key_formats : {key: format(, ...)}
            dict of formats to pre-format the related values
            before insertion into the template.
        """
        super(FormatTemplate, self).__init__()
        self.separator = separator
        self.remove_unused = remove_unused
        self._used_keys = []

        self.key_formats = defaultdict(lambda: ['{:s}'])
        if key_formats:
            self.key_formats.update(key_formats)

    def format(self, format_string, **kwargs):
        """Perform the formatting

        Parameters
        ----------
        format_string : str
            The string to be formatted

        kwargs : dict
            The key/value pairs to insert into the string

        Returns
        -------
        formatted : str
            The formatted string.
        """
        self._used_keys = []

        # Preformat the values
        formatted_kwargs = {}
        for key, value in kwargs.items():
            if value is not None:
                for key_format in self.key_formats[key]:

                    # Get the formatting type character. Indices are:
                    #  0: The first replacement field. There should only be one.
                    #  2: Get the format spec.
                    #  -1: Get the last character representing the type.
                    format_type = list(self.parse(key_format))[0][2][-1]

                    try:
                        value = key_format.format(CONVERSION[format_type](value))
                    except ValueError:
                        pass
                    else:
                        break
                else:
                    raise RuntimeError(
                        'No suitable formatting for {key}: {value} found. Given formatting options:'
                        '\n\t{formats}'.format(
                            key=key, value=value, formats=self.key_formats[key]
                        )
                    )
            formatted_kwargs[key] = value
        result = super(FormatTemplate, self).format(
            format_string, **formatted_kwargs
        )

        # Get any unused arguments and simply do the appending
        unused_keys = set(formatted_kwargs).difference(self._used_keys)
        unused_values = [
            formatted_kwargs[unused]
            for unused in unused_keys
            if formatted_kwargs[unused] is not None
        ]
        result_parts = [result] + unused_values
        result = self.separator.join(result_parts)

        return result

    # Make the instance callable
    __call__ = format

    def get_value(self, key, args, kwargs):
        """Return a given field value

        Parameters
        ----------
        key : str
            The key to retrieve.

        args : [arg(, ...)]
            Positional arguments passed.
            This is ignored.

        kwargs : {k:v(, ...)}
            The key/value pairs passed in.

        Returns
        -------
        obj
            The value from the kwargs.
            If not found, the string '{key}' is returned.
        """
        if self.remove_unused:
            default = ''
        else:
            default = '{' + key + '}'
        value = kwargs.get(key, default)
        self._used_keys.append(key)
        if value is None:
            value = ''

        return value
