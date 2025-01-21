# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
STScI edits to astropy FitsDiff.
"""

import fnmatch
import textwrap

from astropy import __version__
from astropy.utils.diff import diff_values

from astropy.io.fits.card import BLANK_CARD, Card

from astropy.io.fits.hdu.hdulist import HDUList
from astropy.io.fits.hdu.table import _TableLikeHDU

import numpy as np
from astropy.io.fits.diff import FITSDiff, HDUDiff, HeaderDiff, TableDataDiff, ImageDataDiff, RawDataDiff


__all__ = [
    "STFITSDiff",
    "STHDUDiff",
    "STHeaderDiff",
    "ImageDataDiff",
    "RawDataDiff",
    "TableDataDiff",
]


class STFITSDiff(FITSDiff):
    """FITSDiff class from astropy with the STScI edits.

    Diff two FITS files by filename, or two `HDUList` objects

    `FITSDiff` objects have the following diff attributes:

    - ``diff_hdu_count``: If the FITS files being compared have different
      numbers of HDUs, this contains a 2-tuple of the number of HDUs in each
      file.

    - ``diff_hdus``: If any HDUs with the same index are different, this
      contains a list of 2-tuples of the HDU index and the `HDUDiff` object
      representing the differences between the two HDUs.
    """

    def __init__(
        self,
        a,
        b,
        ignore_hdus=[],
        ignore_keywords=[],
        ignore_comments=[],
        ignore_fields=[],
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
        report_pixel_loc_diffs=False,
        extension_tolerances=None,
    ):
        """
        Parameters
        ----------
        a : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object.

        b : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object to
            compare to the first file.

        ignore_hdus : sequence, optional
            HDU names to ignore when comparing two FITS files or HDU lists; the
            presence of these HDUs and their contents are ignored.  Wildcard
            strings may also be included in the list.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        ignore_fields : sequence, optional
            The (case-insensitive) names of any table columns to ignore if any
            table data is to be compared.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found

        extension_tolerances : dict, optional
            Provide a different relative and absolute tolerance for the given extensions, e.g.
            extension_tolerances = {'sci': {'rtol': 1e-3, 'atol': 1e-2},
                                    'err': {'rtol': 1e-1, 'atol': 1e-2},
                                    'VAR_RNOISE': {'rtol': 2, 'atol': 1}}

        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.extension_tolerances = {}
        if extension_tolerances is not None:
            # Make sure the given dict keys are all upper case
            for key in extension_tolerances:
                self.extension_tolerances[key.upper()] = extension_tolerances[key]
            # Make sure the other extensions get a default relative and absolute tolerance
            if 'DEFAULT' not in [key.upper() for key in extension_tolerances]:
                self.extension_tolerances['DEFAULT'] = {'rtol': rtol, 'atol': atol}
        super().__init__(a, b,
                         ignore_hdus=ignore_hdus,
                         ignore_keywords=ignore_keywords,
                         ignore_comments=ignore_comments,
                         ignore_fields=ignore_fields,
                         numdiffs=numdiffs,
                         rtol=rtol,
                         atol=atol,
                         ignore_blanks=ignore_blanks,
                         ignore_blank_cards=ignore_blank_cards)

    def _diff(self):
        if len(self.a) != len(self.b):
            self.diff_hdu_count = (len(self.a), len(self.b))

        # Record filenames for use later in _report
        self.filenamea = self.a.filename()
        if not self.filenamea:
            self.filenamea = f"<{self.a.__class__.__name__} object at {id(self.a):#x}>"

        self.filenameb = self.b.filename()
        if not self.filenameb:
            self.filenameb = f"<{self.b.__class__.__name__} object at {id(self.b):#x}>"

        if self.ignore_hdus:
            self.a = HDUList([h for h in self.a if h.name not in self.ignore_hdus])
            self.b = HDUList([h for h in self.b if h.name not in self.ignore_hdus])
        if self.ignore_hdu_patterns:
            a_names = [hdu.name for hdu in self.a]
            b_names = [hdu.name for hdu in self.b]
            for pattern in self.ignore_hdu_patterns:
                a_ignored = fnmatch.filter(a_names, pattern)
                self.a = HDUList([h for h in self.a if h.name not in a_ignored])
                b_ignored = fnmatch.filter(b_names, pattern)
                self.b = HDUList([h for h in self.b if h.name not in b_ignored])

        # For now, just compare the extensions one by one in order.
        # Might allow some more sophisticated types of diffing later.

        # TODO: Somehow or another simplify the passing around of diff
        # options--this will become important as the number of options grows
        for idx in range(min(len(self.a), len(self.b))):
            if self.extension_tolerances:
                if self.a[idx].name in self.extension_tolerances:
                    self.rtol = self.extension_tolerances[self.a[idx].name]["rtol"]
                    self.atol = self.extension_tolerances[self.a[idx].name]["atol"]
                else:
                    self.rtol = self.extension_tolerances["DEFAULT"]["rtol"]
                    self.atol = self.extension_tolerances["DEFAULT"]["atol"]
            hdu_diff = STHDUDiff.fromdiff(self, self.a[idx], self.b[idx])

            if not hdu_diff.identical:
                if (
                    self.a[idx].name == self.b[idx].name
                    and self.a[idx].ver == self.b[idx].ver
                ):
                    self.diff_hdus.append(
                        (idx, hdu_diff, self.a[idx].name, self.a[idx].ver)
                    )
                else:
                    self.diff_hdus.append((idx, hdu_diff, "", self.a[idx].ver))

    def _report(self):
        wrapper = textwrap.TextWrapper(initial_indent="  ", subsequent_indent="  ")

        self._fileobj.write("\n")
        self._writeln(f" fitsdiff: {__version__}")
        self._writeln(f" a: {self.filenamea}\n b: {self.filenameb}")

        if self.ignore_hdus:
            ignore_hdus = " ".join(sorted(self.ignore_hdus))
            self._writeln(" HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdus))

        if self.ignore_hdu_patterns:
            ignore_hdu_patterns = " ".join(sorted(self.ignore_hdu_patterns))
            self._writeln(
                " HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdu_patterns)
            )

        if self.ignore_keywords:
            ignore_keywords = " ".join(sorted(self.ignore_keywords))
            self._writeln(
                " Keyword(s) not to be compared:\n" + wrapper.fill(ignore_keywords)
            )

        if self.ignore_comments:
            ignore_comments = " ".join(sorted(self.ignore_comments))
            self._writeln(
                " Keyword(s) whose comments are not to be compared:\n"
                + wrapper.fill(ignore_comments)
            )

        if self.ignore_fields:
            ignore_fields = " ".join(sorted(self.ignore_fields))
            self._writeln(
                " Table column(s) not to be compared:\n" + wrapper.fill(ignore_fields)
            )

        self._writeln(
            f" Maximum number of different data values to be reported: {self.numdiffs}"
        )

        if not self.extension_tolerances:
            self._writeln(f"\n Relative tolerance: {self.rtol}, Absolute tolerance: {self.atol}")

        if self.diff_hdu_count:
            self._fileobj.write("\n")
            self._writeln("Files contain different numbers of HDUs:")
            self._writeln(f" a: {self.diff_hdu_count[0]}")
            self._writeln(f" b: {self.diff_hdu_count[1]}")

            if not self.diff_hdus:
                self._writeln("No differences found between common HDUs.")
                return
        elif not self.diff_hdus:
            self._fileobj.write("\n")
            self._writeln("No differences found.")
            return

        for idx, hdu_diff, extname, extver in self.diff_hdus:
            # print out the extension heading
            if not self.extension_tolerances:
                if idx == 0:
                    self._fileobj.write("\n")
                    self._writeln("Primary HDU:")
                else:
                    self._fileobj.write("\n")
                    if extname:
                        self._writeln(f"Extension HDU {idx} ({extname}, {extver}):")
                    else:
                        self._writeln(f"Extension HDU {idx}:")
            else:
                self._writeln(f"Extension HDU {idx} ({extname}, {extver}):")
                if extname in self.extension_tolerances:
                    rtol = self.extension_tolerances[extname]["rtol"]
                    atol = self.extension_tolerances[extname]["atol"]
                else:
                    rtol = self.extension_tolerances["DEFAULT"]["rtol"]
                    atol = self.extension_tolerances["DEFAULT"]["atol"]
                self._writeln(
                    f"\n Relative tolerance: {rtol}, Absolute tolerance: {atol}"
                )
            hdu_diff.report(self._fileobj, indent=self._indent + 1)


class STHDUDiff(HDUDiff):
    """
    HDUDiff class from astropy with the STScI edits.

    Diff two HDU objects, including their headers and their data (but only if
    both HDUs contain the same type of data (image, table, or unknown).

    `HDUDiff` objects have the following diff attributes:

    - ``diff_extnames``: If the two HDUs have different EXTNAME values, this
      contains a 2-tuple of the different extension names.

    - ``diff_extvers``: If the two HDUS have different EXTVER values, this
      contains a 2-tuple of the different extension versions.

    - ``diff_extlevels``: If the two HDUs have different EXTLEVEL values, this
      contains a 2-tuple of the different extension levels.

    - ``diff_extension_types``: If the two HDUs have different XTENSION values,
      this contains a 2-tuple of the different extension types.

    - ``diff_headers``: Contains a `HeaderDiff` object for the headers of the
      two HDUs. This will always contain an object--it may be determined
      whether the headers are different through ``diff_headers.identical``.

    - ``diff_data``: Contains either a `ImageDataDiff`, `TableDataDiff`, or
      `RawDataDiff` as appropriate for the data in the HDUs, and only if the
      two HDUs have non-empty data of the same type (`RawDataDiff` is used for
      HDUs containing non-empty data of an indeterminate type).
    """

    def __init__(
        self,
        a,
        b,
        ignore_keywords=[],
        ignore_comments=[],
        ignore_fields=[],
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
        report_pixel_loc_diffs=False,
    ):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        ignore_fields : sequence, optional
            The (case-insensitive) names of any table columns to ignore if any
            table data is to be compared.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.nans = []
        self.percentages = {}
        self.stats = {}
        super().__init__(a, b,
                         ignore_keywords=ignore_keywords,
                         ignore_comments=ignore_comments,
                         ignore_fields=ignore_fields,
                         numdiffs=numdiffs,
                         rtol=rtol,
                         atol=atol,
                         ignore_blanks=ignore_blanks,
                         ignore_blank_cards=ignore_blank_cards)

    def _diff(self):
        if self.a.name != self.b.name:
            self.diff_extnames = (self.a.name, self.b.name)

        if self.a.ver != self.b.ver:
            self.diff_extvers = (self.a.ver, self.b.ver)

        if self.a.level != self.b.level:
            self.diff_extlevels = (self.a.level, self.b.level)

        if self.a.header.get("XTENSION") != self.b.header.get("XTENSION"):
            self.diff_extension_types = (
                self.a.header.get("XTENSION"),
                self.b.header.get("XTENSION"),
            )

        self.diff_headers = STHeaderDiff.fromdiff(
            self, self.a.header.copy(), self.b.header.copy()
        )

        def get_percentage(values, threshold):
            # Calculate percentage
            n_total = values.size
            n = np.logical_or(values > threshold, values < -threshold).nonzero()[0].size
            return np.round((n / n_total) * 100, decimals=2)

        def get_quick_report(a, b):
            # match nans for all arrays and remove them for logical comparison
            nan_idx = (np.isnan(a) | np.isnan(b))
            anonan = a[~nan_idx]
            bnonan = b[~nan_idx]
            # Get the number of NaN in each array
            nans = [np.isnan(a).size, np.isnan(b).size]
            # Calculate stats
            values = np.abs(np.abs(anonan) - np.abs(bnonan))
            stats = {'mean_value_in_a': np.mean(anonan),
                     'mean_value_in_b': np.mean(bnonan),
                     'max_diff': max(values),
                     'min_diff': min(values),
                     'mean_diff': np.mean(values),
                     'std_dev_diff': np.std(values)}
            # Calculate difference percentages
            percentages = {}
            thresholds = [0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0]
            for threshold in thresholds:
                percentage = get_percentage(values, threshold)
                percentages[threshold] = percentage
            return nans, percentages, stats

        if self.a.data is None or self.b.data is None:
            # TODO: Perhaps have some means of marking this case
            pass
        elif self.a.is_image and self.b.is_image:
            self.diff_data = ImageDataDiff.fromdiff(self, self.a.data, self.b.data)
            self.nans, self.percentages, self.stats = get_quick_report(self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif isinstance(self.a, _TableLikeHDU) and isinstance(self.b, _TableLikeHDU):
            # TODO: Replace this if/when _BaseHDU grows a .is_table property
            self.diff_data = TableDataDiff.fromdiff(self, self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif not self.diff_extension_types:
            # Don't diff the data for unequal extension types that are not
            # recognized image or table types
            self.diff_data = RawDataDiff.fromdiff(self, self.a.data, self.b.data)
            self.nans, self.percentages, self.stats = get_quick_report(self.a.data, self.b.data)

            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None

    def _report(self):
        if self.identical:
            self._writeln(" No differences found.")
        if self.diff_extension_types:
            self._writeln(
                " Extension types differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extension_types
                )
            )
        if self.diff_extnames:
            self._writeln(
                " Extension names differ:\n  a: {}\n  b: {}".format(*self.diff_extnames)
            )
        if self.diff_extvers:
            self._writeln(
                " Extension versions differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extvers
                )
            )

        if self.diff_extlevels:
            self._writeln(
                " Extension levels differ:\n  a: {}\n  b: {}".format(
                    *self.diff_extlevels
                )
            )

        if not self.diff_headers.identical:
            self._fileobj.write("\n")
            self._writeln(" Headers contain differences:")
            self.diff_headers.report(self._fileobj, indent=self._indent + 1)

        def report_data_diff():
            self._writeln("  NaN in arrays:")
            if len(self.nans) > 0:
                self._writeln(f"  a: {self.nans[0]}")
                self._writeln(f"  b: {self.nans[1]}")
            else:
                self._writeln("  a and b: 0")
            # Calculate difference percentages
            self._writeln(" Difference of a from b:")
            for key, val in self.percentages.items():
                self._writeln(f"  {key:>6} ..... {val:<5}%")
            self._writeln(" Stats:")
            for key, val in self.stats.items():
                self._writeln(f"  {key} = {val}")

        if self.diff_data is not None and not self.diff_data.identical:
            self._fileobj.write("\n")
            self._writeln(" Data contains differences:")
            report_data_diff()
            if self.report_pixel_loc_diffs:
                self.diff_data.report(self._fileobj, indent=self._indent + 1)


class STHeaderDiff(HeaderDiff):
    """
    HeaderDiff class from astropy with the STScI edits.

    Diff two `Header` objects.

    `HeaderDiff` objects have the following diff attributes:

    - ``diff_keyword_count``: If the two headers contain a different number of
      keywords, this contains a 2-tuple of the keyword count for each header.

    - ``diff_keywords``: If either header contains one or more keywords that
      don't appear at all in the other header, this contains a 2-tuple
      consisting of a list of the keywords only appearing in header a, and a
      list of the keywords only appearing in header b.

    - ``diff_duplicate_keywords``: If a keyword appears in both headers at
      least once, but contains a different number of duplicates (for example, a
      different number of HISTORY cards in each header), an item is added to
      this dict with the keyword as the key, and a 2-tuple of the different
      counts of that keyword as the value.  For example::

          {'HISTORY': (20, 19)}

      means that header a contains 20 HISTORY cards, while header b contains
      only 19 HISTORY cards.

    - ``diff_keyword_values``: If any of the common keyword between the two
      headers have different values, they appear in this dict.  It has a
      structure similar to ``diff_duplicate_keywords``, with the keyword as the
      key, and a 2-tuple of the different values as the value.  For example::

          {'NAXIS': (2, 3)}

      means that the NAXIS keyword has a value of 2 in header a, and a value of
      3 in header b.  This excludes any keywords matched by the
      ``ignore_keywords`` list.

    - ``diff_keyword_comments``: Like ``diff_keyword_values``, but contains
      differences between keyword comments.

    `HeaderDiff` objects also have a ``common_keywords`` attribute that lists
    all keywords that appear in both headers.
    """

    def __init__(
        self,
        a,
        b,
        ignore_keywords=[],
        ignore_comments=[],
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
    ):
        """
        Parameters
        ----------
        a : `~astropy.io.fits.Header` or string or bytes
            A header.

        b : `~astropy.io.fits.Header` or string or bytes
            A header to compare to the first header.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers; the presence
            of these keywords and their values are ignored.  Wildcard strings
            may also be included in the list.

        ignore_comments : sequence, optional
            A list of header keywords whose comments should be ignored in the
            comparison.  May contain wildcard strings as with ignore_keywords.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        rtol : float, optional
            The relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns
            (default: 0.0). Values which satisfy the expression

            .. math::

                \\left| a - b \\right| > \\text{atol} + \\text{rtol} \\cdot \\left| b \\right|

            are considered to be different.
            The underlying function used for comparison is `numpy.allclose`.

            .. versionadded:: 2.0

        atol : float, optional
            The allowed absolute difference. See also ``rtol`` parameter.

            .. versionadded:: 2.0

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored
            (default: True).

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain
            whitespace (default: True).
        """
        super().__init__(a, b,
                         ignore_keywords=ignore_keywords,
                         ignore_comments=ignore_comments,
                         rtol=rtol,
                         atol=atol,
                         ignore_blanks=ignore_blanks,
                         ignore_blank_cards=ignore_blank_cards)

    def _diff(self):
        if self.ignore_blank_cards:
            cardsa = [c for c in self.a.cards if str(c) != BLANK_CARD]
            cardsb = [c for c in self.b.cards if str(c) != BLANK_CARD]
        else:
            cardsa = list(self.a.cards)
            cardsb = list(self.b.cards)

        # build dictionaries of keyword values and comments
        def get_header_values_comments(cards):
            values = {}
            comments = {}
            for card in cards:
                value = card.value
                if self.ignore_blanks and isinstance(value, str):
                    value = value.rstrip()
                values.setdefault(card.keyword, []).append(value)
                comments.setdefault(card.keyword, []).append(card.comment)
            return values, comments

        valuesa, commentsa = get_header_values_comments(cardsa)
        valuesb, commentsb = get_header_values_comments(cardsb)

        keywordsa = set(valuesa)
        keywordsb = set(valuesb)

        self.common_keywords = sorted(keywordsa.intersection(keywordsb))
        if len(cardsa) != len(cardsb):
            self.diff_keyword_count = (len(cardsa), len(cardsb))

        # Any other diff attributes should exclude ignored keywords
        keywordsa = keywordsa.difference(self.ignore_keywords)
        keywordsb = keywordsb.difference(self.ignore_keywords)
        if self.ignore_keyword_patterns:
            for pattern in self.ignore_keyword_patterns:
                keywordsa = keywordsa.difference(fnmatch.filter(keywordsa, pattern))
                keywordsb = keywordsb.difference(fnmatch.filter(keywordsb, pattern))

        if "*" in self.ignore_keywords:
            # Any other differences between keywords are to be ignored
            return

        left_only_keywords = sorted(keywordsa.difference(keywordsb))
        right_only_keywords = sorted(keywordsb.difference(keywordsa))

        if left_only_keywords or right_only_keywords:
            self.diff_keywords = (left_only_keywords, right_only_keywords)

        # Compare count of each common keyword
        for keyword in self.common_keywords:
            if keyword in self.ignore_keywords:
                continue
            if self.ignore_keyword_patterns:
                skip = False
                for pattern in self.ignore_keyword_patterns:
                    if fnmatch.fnmatch(keyword, pattern):
                        skip = True
                        break
                if skip:
                    continue

            counta = len(valuesa[keyword])
            countb = len(valuesb[keyword])
            if counta != countb:
                self.diff_duplicate_keywords[keyword] = (counta, countb)

            # Compare keywords' values and comments
            for a in valuesa[keyword]:
                if a not in valuesb[keyword]:
                    continue
                bidx = valuesb[keyword].index(a)
                b = valuesb[keyword][bidx]

                if diff_values(a, b, rtol=self.rtol, atol=self.atol):
                    self.diff_keyword_values[keyword].append((a, b))
                else:
                    # If there are duplicate keywords we need to be able to
                    # index each duplicate; if the values of a duplicate
                    # are identical use None here
                    self.diff_keyword_values[keyword].append(None)

            if not any(self.diff_keyword_values[keyword]):
                # No differences found; delete the array of Nones
                del self.diff_keyword_values[keyword]

            if "*" in self.ignore_comments or keyword in self.ignore_comments:
                continue
            if self.ignore_comment_patterns:
                skip = False
                for pattern in self.ignore_comment_patterns:
                    if fnmatch.fnmatch(keyword, pattern):
                        skip = True
                        break
                if skip:
                    continue

            for a in commentsa[keyword]:
                if a not in commentsb[keyword]:
                    continue
                bidx = commentsb[keyword].index(a)
                b = commentsb[keyword][bidx]
                if diff_values(a, b):
                    self.diff_keyword_comments[keyword].append((a, b))
                else:
                    self.diff_keyword_comments[keyword].append(None)

            if not any(self.diff_keyword_comments[keyword]):
                del self.diff_keyword_comments[keyword]


