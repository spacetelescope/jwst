"""
A module that provides functions for manipulating bitmasks and data quality (DQ) arrays.

:Authors: Mihai Cara (contact: help@stsci.edu)

:License: `<http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE>`_

"""

import numpy as np

__version__ = '0.1.0'
__vdate__ = '29-March-2015'
__author__ = 'Mihai Cara'


def interpret_bits_value(val):
    """
    Converts input bits value from string to a single integer value or None.
    If a comma- or '+'-separated set of values are provided, they are summed.

    .. note::
        In order to flip the bits of the final result (after summation),
        for input of `str` type, prepend '~' to the input string. '~' must
        be prepended to the *entire string* and not to each bit flag!

    Parameters
    ----------
    val : int, str, None
        An integer bit mask or flag, `None`, or a comma- or '+'-separated
        string list of integer bit values. If `val` is a `str` and if
        it is prepended with '~', then the output bit mask will have its
        bits flipped (compared to simple sum of input val).

    Returns
    -------
    bitmask : int or None
        Returns and integer bit mask formed from the input bit value
        or `None` if input `val` parameter is `None` or an empty string.
        If input string value was prepended with '~', then returned
        value will have its bits flipped (inverse mask).

    Examples
    --------
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags(28) )
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('4,8,16') )
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('~4,8,16') )
        '1111111111100011'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('~(4+8+16)') )
        '1111111111100011'

    """
    if isinstance(val, int) or val is None:
        return val

    else:
        val = str(val).strip()

        if val.startswith('~'):
            flip_bits = True
            val = val[1:].lstrip()
        else:
            flip_bits = False

        if val.startswith('('):
            if val.endswith(')'):
                val = val[1:-1].strip()
            else:
                raise ValueError('Unbalanced parantheses or incorrect syntax.')

        if ',' in val:
            valspl = val.split(',')
            bitmask = 0
            for v in valspl:
                bitmask += int(v)

        elif '+' in val:
            valspl = val.split('+')
            bitmask = 0
            for v in valspl:
                bitmask += int(v)

        elif val.upper() in ['', 'NONE', 'INDEF']:
            return None

        else:
            bitmask = int(val)

        if flip_bits:
            bitmask = ~bitmask

    return bitmask


def bitmask2mask(bitmask, ignore_bits, good_mask_value=1, dtype=np.uint8):
    """
    bitmask2mask(bitmask, ignore_bits, good_pix_value=1, dtype=numpy.uint8)
    Interprets an array of bit flags and converts it to a "binary" mask array.
    This function is particularly useful to convert data quality arrays to
    binary masks.

    Parameters
    ----------
    bitmask : numpy.ndarray
        An array of bit flags. Values different from zero are interpreted as
        "bad" values and values equal to zero are considered as "good" values.
        However, see `ignore_bits` parameter on how to ignore some bits
        in the `bitmask` array.

    ignore_bits : int, str, None
        An integer bit mask, `None`, or a comma- or '+'-separated
        string list of integer bit values that indicate what bits in the
        input `bitmask` should be *ignored* (i.e., zeroed). If `ignore_bits`
        is a `str` and if it is prepended with '~', then the meaning
        of `ignore_bits` parameters will be reversed: now it will be
        interpreted as a list of bits to be *used* (or *not ignored*) when
        deciding what elements of the input `bitmask` array are "bad".

        The `ignore_bits` parameter is the integer sum of all of the bit
        values from the input `bitmask` array that should be considered
        "good" when creating the output binary mask. For example, if
        values in the `bitmask` array can be combinations
        of 1, 2, 4, and 8 flags and one wants to consider that
        values having *only* bit flags 2 and/or 4 as being "good",
        then `ignore_bits` should be set to 2+4=6. Then a `bitmask` element
        having values 2,4, or 6 will be considered "good", while an
        element with a value, e.g., 1+2=3, 4+8=12, etc. will be interpreted
        as "bad".

        Alternatively, one can enter a comma- or '+'-separated list
        of integer bit flags that should be added to obtain the
        final "good" bits. For example, both ``4,8`` and ``4+8``
        are equivalent to setting `ignore_bits` to 12.

        See :py:func:`interpret_bits_value` for examples.

        | Setting `ignore_bits` to `None` effectively will interpret
          all `bitmask` elements as "good" regardless of their value.

        | Setting `ignore_bits` to 0 effectively will assume that all
          non-zero elements in the input `bitmask` array are to be
          interpreted as "bad".

        | In order to reverse the meaning of the `ignore_bits`
          parameter from indicating bits in the values of `bitmask`
          elements that should be ignored when deciding which elements
          are "good" (these are the elements that are zero after ignoring
          `ignore_bits`), to indicating the bits should be used
          exclusively in deciding whether a `bitmask` element is "good",
          prepend '~' to the string value. For example, in order to use
          **only** (or **exclusively**) flags 4 and 8 (2nd and 3rd bits)
          in the values of the input `bitmask` array when deciding whether
          or not that element is "good", set `ignore_bits` to ``~4+8``,
          or ``~4,8 To obtain the same effect with an `int` input value
          (except for 0), enter -(4+8+1)=-9. Following this convention,
          a `ignore_bits` string value of ``'~0'`` would be equivalent to
          setting ``ignore_bits=None``.

    good_mask_value : int, bool (Default = 1)
        This parameter is used to derive the values that will be assigned to
        the elements in the output `mask` array that correspond to the "good"
        flags (that are 0 after zeroing bits specified by `ignore_bits`)
        in the input `bitmask` array. When `good_mask_value` is non-zero or
        `True` then values in the output mask array corresponding to "good"
        bit flags in `bitmask` will be 1 (or `True` if `dtype` is `bool`) and
        values of corresponding to "bad" flags will be 0. When
        `good_mask_value` is zero or `False` then values in the output mask
        array corresponding to "good" bit flags in `bitmask` will be 0
        (or `False` if `dtype` is `bool`) and values of corresponding
        to "bad" flags will be 1.

    dtype : data-type (Default = numpy.uint8)
        The desired data-type for the output binary mask array.

    Returns
    -------
    mask : numpy.ndarray
        Returns an array whose elements can have two possible values,
        e.g., 1 or 0 (or `True` or `False` if `dtype` is `bool`) according to
        values of to the input `bitmask` elements, `ignore_bits` parameter,
        and the `good_mask_value` parameter.

    Examples
    --------
        >>> from stsci.tools import bitmask
        >>> dqbits = np.asarray([[0,0,1,2,0,8,12,0],[10,4,0,0,0,16,6,0]])
        >>> bitmask.bitmask2mask(dqbits, ignore_bits=0, dtype=int)
        array([[1, 1, 0, 0, 1, 0, 0, 1],
               [0, 0, 1, 1, 1, 0, 0, 1]])
        >>> bitmask.bitmask2mask(dqbits, ignore_bits=0, dtype=bool)
        array([[ True,  True, False, False,  True, False, False,  True],
               [False, False,  True,  True,  True, False, False,  True]], dtype=bool)
        >>> bitmask.bitmask2mask(dqbits, ignore_bits=6, good_pix_value=0, dtype=int)
        array([[0, 0, 1, 0, 0, 1, 1, 0],
               [1, 0, 0, 0, 0, 1, 0, 0]])
        >>> bitmask.bitmask2mask(dqbits, ignore_bits=~6, good_pix_value=0, dtype=int)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitmask2mask(dqbits, ignore_bits='~(2+4)', good_pix_value=0, dtype=int)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])

    """

    ignore_bits = interpret_bits_value(ignore_bits)

    if good_mask_value:
        mask = np.ones_like(bitmask, dtype=dtype)
        if ignore_bits is None:
            return mask
        bad_mask_value = 0

    else:
        mask = np.zeros_like(bitmask, dtype=dtype)
        if ignore_bits is None:
            return mask
        bad_mask_value = 1

    mask[np.bitwise_and(bitmask, ~ignore_bits) > 0] = bad_mask_value

    return mask
