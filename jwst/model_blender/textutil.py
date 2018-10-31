"""Text output-related utilities."""


import textwrap


def wrap(text, width, *args, **kwargs):
    """
    Like :func:`textwrap.wrap` but preserves existing newlines which
    :func:`textwrap.wrap` does not otherwise handle well.

    See Also
    --------
    :func:`textwrap.wrap`
    """

    return sum([textwrap.wrap(line, width, *args, **kwargs)
                if line else [''] for line in text.splitlines()], [])


def textbox(text, width=78, boxchar='#', indent=0):
    """
    Outputs line-wrapped text wrapped in a box drawn with a repeated (usually
    ASCII) character.

    For example:

        >>> print(textbox('Text to wrap', width=16))
        ################
        #              #
        # Text to wrap #
        #              #
        ################

    Parameters
    ----------
    text : string
        The text to wrap

    width : int
        The width of the entire box, including the perimeter and
        the indentation space.  Because the
        wrapped text is padded with an additional column of whitespace on each
        side, the minimum width is 5--any width less than that is
        is automatically increased to 5 (default: 78)

    boxchar : string
        (No pun intended.) The character to draw the box with.  May also
        be a string of multiple characters (default: '#')

    indent : int
        Amount of space by which the box should be indented. (default: 0)

    """

    min_width = len(boxchar) * 2 + 3
    width = max(width-indent, min_width)
    indentspace = indent * ' '

    wrap_width = width - min_width + 1

    q, r = divmod(width, len(boxchar))
    # The top/bottom border
    top_border = indentspace + boxchar * q + boxchar[:r]
    top_padding  = indentspace + boxchar + ' ' * (width - len(boxchar) * 2) + boxchar

    lines = ['%s%s %s %s' % (indentspace, boxchar, line.ljust(wrap_width),
                             boxchar)
             for line in wrap(text, wrap_width)]
    top = [top_border, top_padding]
    bottom = [top_padding, top_border]
    return '\n'.join(top + lines + bottom)
