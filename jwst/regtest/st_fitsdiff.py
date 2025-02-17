# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
STScI edits to astropy FitsDiff.
"""

import fnmatch
import operator
import textwrap
from itertools import islice

from astropy import __version__
from astropy.utils.diff import diff_values, report_diff_values, where_not_allclose

from astropy.io.fits.card import BLANK_CARD, Card

from astropy.io.fits.hdu.hdulist import HDUList
from astropy.io.fits.hdu.table import _TableLikeHDU

import numpy as np
from astropy.io.fits.diff import (FITSDiff, HDUDiff, HeaderDiff, TableDataDiff, ImageDataDiff,
                                  report_diff_keyword_attr, _get_differences, _COL_ATTRS)


__all__ = [
    "STFITSDiff",
    "STHDUDiff",
    "STHeaderDiff",
    "STImageDataDiff",
    "STRawDataDiff",
    "STTableDataDiff",
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
                                    'VAR_RNOISE': {'rtol': 2, 'atol': 1},
                                    'default': {'rtol': 1e-5, 'atol': 1e-7}}
            The key for only providing main header tolerances is 'primary'.
            The key 'headers' can be used to provide a special tolerance for all
            extension headers in the HDU list.
            It does not matter if the keys in the dictionary are upper or lower case.
            The key 'default' is optional, i.e. if it is not provided then the default values will
            be used, otherwise the default value will be the one in the dictionary.

        """
        self.diff_dimensions = ()
        self.diff_extnames = ()
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.header_tolerances = {}
        self.extension_tolerances = {}
        if extension_tolerances is not None:
            # Make sure the given dict keys are all upper case
            for key in extension_tolerances:
                # Make sure both tolerances exist in key
                tols = extension_tolerances[key]
                if 'rtol' not in tols:
                    tols['rtol'] = rtol
                if 'atol' not in tols:
                    tols['atol'] = atol
                if isinstance(key, str):
                    self.extension_tolerances[key.upper()] = tols
                else:
                    self.extension_tolerances[key] = tols
            # Make sure the other extensions get a default relative and absolute tolerance
            if 'DEFAULT' not in [key.upper() for key in extension_tolerances if isinstance(key, str)]:
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

        # Check if the files contain the same extension names and find the intersection
        ext_namesa = [ext.name for ext in self.a]
        ext_namesb = [ext.name for ext in self.b]
        ext_not_in_both = set(ext_namesa).symmetric_difference(set(ext_namesb))
        ext_intersection = set.intersection(set(ext_namesa), set(ext_namesb))
        if self.diff_hdu_count:
            self.diff_extnames = (ext_namesa, ext_namesb, ext_intersection, ext_not_in_both)

        # Set the tolerance for the headers
        if "HEADERS" in [key for key in self.extension_tolerances]:
            self.header_tolerances["rtol"] = self.extension_tolerances["HEADERS"]["rtol"]
            self.header_tolerances["atol"] = self.extension_tolerances["HEADERS"]["atol"]

        # Make sure to compare the same extensions
        for idxa, extname in enumerate(ext_namesa):
            extver = self.a[idxa].ver
            same_version = False
            if extname not in ext_intersection:
                continue
            else:
                # make sure to compare same version of extension
                idx_extname_b = [i for i, x in enumerate(ext_namesb) if x == extname]
                for idxb in idx_extname_b:
                    if extver == self.b[idxb].ver:
                        same_version = True
                        break
                if same_version:
                    if self.extension_tolerances:
                        if idxa in self.extension_tolerances:
                            if idxa in self.extension_tolerances:
                                self.rtol = self.extension_tolerances[idxa]["rtol"]
                                self.atol = self.extension_tolerances[idxa]["atol"]
                            elif extname in self.extension_tolerances:
                                self.rtol = self.extension_tolerances[extname]["rtol"]
                                self.atol = self.extension_tolerances[extname]["atol"]
                            else:
                                self.rtol = self.extension_tolerances["DEFAULT"]["rtol"]
                                self.atol = self.extension_tolerances["DEFAULT"]["atol"]
                    hdu_diff = STHDUDiff.fromdiff(self, self.a[idxa], self.b[idxb])
                    if not hdu_diff.identical:
                        self.diff_hdus.append((idxa, hdu_diff, extname, extver))

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
            self._writeln(f" a: {self.diff_hdu_count[0]}, {self.diff_extnames[0]}")
            self._writeln(f" b: {self.diff_hdu_count[1]}, {self.diff_extnames[1]}")
            self._writeln(f" Common HDUs: {[en for en in self.diff_extnames[2]]}")
            self._writeln(f" Missing HDUs: {[en for en in self.diff_extnames[3]]}")

            if not self.diff_hdus:
                self._fileobj.write("\n")
                self._writeln("No differences found between common HDUs.")
                return

        elif not self.diff_hdus:
            self._fileobj.write("\n")
            self._writeln("No differences found.")
            return

        self._fileobj.write("\n")
        for idx, hdu_diff, extname, extver in self.diff_hdus:
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
                self._fileobj.write("\n")
                self._writeln(f"Extension HDU {idx} ({extname}, {extver}):")
                if idx in self.extension_tolerances:
                    rtol = self.extension_tolerances[idx]["rtol"]
                    atol = self.extension_tolerances[idx]["atol"]
                elif extname in self.extension_tolerances:
                    rtol = self.extension_tolerances[extname]["rtol"]
                    atol = self.extension_tolerances[extname]["atol"]
                else:
                    rtol = self.extension_tolerances["DEFAULT"]["rtol"]
                    atol = self.extension_tolerances["DEFAULT"]["atol"]
                self._writeln(
                    f"\n  Relative tolerance: {rtol:.1e}, Absolute tolerance: {atol:.1e}"
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
        header_tolerances={}
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

        header_tolerances : dict, optional
            Dictionary with the relative and absolute tolerances for all headers.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.header_tolerances = header_tolerances
        self.nans, self.percentages, self.stats = None, None, None
        self.diff_dimensions = ()
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

        # If specific header tolerances were given set them temporarily
        if self.header_tolerances:
            rtol, atol = self.rtol, self.atol
            self.rtol, self.atol = self.header_tolerances['rtol'], self.header_tolerances['atol']

        # Get the header differences
        self.diff_headers = STHeaderDiff.fromdiff(
            self, self.a.header.copy(), self.b.header.copy()
        )
        # Reset the object tolerances
        if self.header_tolerances:
            self.rtol, self.atol = rtol, atol

        def get_quick_report(a, b):
            # Get the number of NaN in each array and other info
            nans_zero_info = [a[np.isnan(a)].size, b[np.isnan(b)].size,
                              a[~np.isnan(a)].size, b[~np.isnan(b)].size,
                              a[a == 0.0].size, b[b == 0.0].size]
            # Match nans for all arrays and remove them for logical comparison
            percentages, stats = {}, {}
            shapea = a.shape
            shapeb = b.shape
            if shapea != shapeb:
                percentages['array_shapes_are_different'] = ''
                stats['no_stats_available'] = ''
                return nans_zero_info, percentages, stats
            nan_idx = (np.isnan(a) | np.isnan(b))
            anonan = a[~nan_idx]
            bnonan = b[~nan_idx]
            values = np.abs(anonan - bnonan)
            if (values == 0.0).all():
                return None, None, None
            # Calculate stats
            stats['mean_value_in_a'] = np.mean(anonan)
            stats['mean_value_in_b'] = np.mean(bnonan)
            # Catch the all NaNs case
            if values.size == 0:
                percentages['NaN'] = 100
                stats['no_stats_available'] = ''
                return nans_zero_info, percentages, stats
            stats['max_abs_diff'] = np.max(values)
            stats['min_abs_diff'] = np.min(values)
            stats['mean_abs_diff'] = np.mean(values)
            stats['std_dev_abs_diff'] = np.std(values)
            nozeros = (values != 0.0) & (bnonan != 0.0)
            relative_values = values[nozeros] / np.abs(bnonan[nozeros])
            # Catch an empty sequence
            if relative_values.size == 0:
                stats['no_rel_stats_available'] = np.nan
            else:
                stats['max_rel_diff'] = np.max(relative_values)
                if 0.0 in values:
                    stats['min_rel_diff'] = 0.0
                else:
                    stats['min_rel_diff'] = np.min(relative_values)
                stats['mean_rel_diff'] = np.mean(relative_values)
                stats['std_dev_rel_diff'] = np.std(relative_values)
            # Calculate difference percentages
            thresholds = [0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0]
            n_total = values.size
            percent_abs_list = []
            for threshold in thresholds:
                n = values[values > threshold + self.atol].size
                percent_abs = np.round((n / n_total) * 100, decimals=2)
                percent_abs_list.append(percent_abs)
            if (percent_abs_list == percent_abs).all():
                percentages['abdiffs_too_large'] = percent_abs
            else:
                # Only include the percentage for 0.0
                percentages['0.0_abs'] = percent_abs
            if relative_values.size > 0:
                percent_rel_list = []
                # Differences are too large values. Calculating percentages on relative numbers.
                for threshold in thresholds:
                    n = relative_values[relative_values > threshold + self.rtol].size
                    percent_rel = np.round((n / n_total) * 100, decimals=2)
                    percentages[str(threshold) + '_rel'] = percent_rel
                    percent_rel_list.append(percent_rel)
                if (percent_rel_list == percent_rel).all():
                    percentages['reldiffs_too_large'] = percent_rel
            return nans_zero_info, percentages, stats

        if self.a.data is None or self.b.data is None:
            # TODO: Perhaps have some means of marking this case
            pass
        elif self.a.is_image and self.b.is_image:
            self.diff_data = STImageDataDiff.fromdiff(self, self.a.data, self.b.data)
            if self.diff_data.diff_total > 0:
                self.nans, self.percentages, self.stats = get_quick_report(self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif isinstance(self.a, _TableLikeHDU) and isinstance(self.b, _TableLikeHDU):
            # TODO: Replace this if/when _BaseHDU grows a .is_table property
            self.diff_data = STTableDataDiff.fromdiff(self, self.a.data, self.b.data)
            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None
        elif not self.diff_extension_types:
            # Don't diff the data for unequal extension types that are not
            # recognized image or table types
            self.diff_data = STRawDataDiff.fromdiff(self, self.a.data, self.b.data)
            if self.diff_data.diff_total > 0:
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

        if self.diff_dimensions:
            dimsa = " x ".join(str(d) for d in reversed(self.diff_dimensions[0]))
            dimsb = " x ".join(str(d) for d in reversed(self.diff_dimensions[1]))
            self._writeln(" Data dimensions differ:")
            self._writeln(f"  a: {dimsa}")
            self._writeln(f"  b: {dimsb}")

        def report_data_diff():
            if len(self.nans) > 0:
                self._writeln(" NaN in arrays:")
                self._writeln(f"  a: {self.nans[0]}")
                self._writeln(f"  b: {self.nans[1]}")
                self._writeln(" No NaN values in arrays:")
                self._writeln(f"  a: {self.nans[2]}")
                self._writeln(f"  b: {self.nans[3]}")
                self._writeln(" Zeros in arrays:")
                self._writeln(f"  a: {self.nans[4]}")
                self._writeln(f"  b: {self.nans[5]}")
            # Calculate difference percentages
            self._writeln(" Difference of a from b:")
            if 'abdiffs_too_large' in self.percentages:
                self._writeln("  * Absolute number differences are too large.")
                self._writeln(f"  {'0.0_abs':>10} ..... {self.percentages['abdiffs_too_large']:<5}%")
                del self.percentages['abdiffs_too_large']
            if 'reldiffs_too_large' in self.percentages:
                self._writeln("  * Relative number differences are too large.")
                self._writeln(f"  {'0.0_rel':>10} ..... {self.percentages['reldiffs_too_large']:<5}%")
            else:
                for key, val in self.percentages.items():
                    self._writeln(f"  {key:>10} ..... {val:<5}%")
            self._writeln(" Stats:")
            for key, val in self.stats.items():
                self._writeln(f"  {key} = {val}")

        if self.diff_data is not None and not self.diff_data.identical:
            self._fileobj.write("\n")
            if [self.nans, self.percentages, self.stats] != [None, None, None]:
                report_data_diff()
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

    def _report(self):
        if self.diff_keyword_count:
            self._writeln(" Headers have different number of cards:")
            self._writeln(f"  a: {self.diff_keyword_count[0]}")
            self._writeln(f"  b: {self.diff_keyword_count[1]}")
        if self.diff_keywords:
            for keyword in self.diff_keywords[0]:
                if keyword in Card._commentary_keywords:
                    val = self.a[keyword][0]
                else:
                    val = self.a[keyword]
                self._writeln(f" Extra keyword {keyword!r:8} in a: {val!r}")
            for keyword in self.diff_keywords[1]:
                if keyword in Card._commentary_keywords:
                    val = self.b[keyword][0]
                else:
                    val = self.b[keyword]
                self._writeln(f" Extra keyword {keyword!r:8} in b: {val!r}")

        if self.diff_duplicate_keywords:
            for keyword, count in sorted(self.diff_duplicate_keywords.items()):
                self._writeln(f" Inconsistent duplicates of keyword {keyword!r:8}:")
                self._writeln(
                    "  Occurs {} times in a, {} times in (b)".format(*count)
                )

        if self.diff_keyword_values or self.diff_keyword_comments:
            for keyword in self.common_keywords:
                report_diff_keyword_attr(
                    self._fileobj,
                    "values",
                    self.diff_keyword_values,
                    keyword,
                    ind=self._indent,
                )
                report_diff_keyword_attr(
                    self._fileobj,
                    "comments",
                    self.diff_keyword_comments,
                    keyword,
                    ind=self._indent,
                )


class STImageDataDiff(ImageDataDiff):
    """
    Diff two image data arrays (really any array from a PRIMARY HDU or an IMAGE
    extension HDU, though the data unit is assumed to be "pixels").

    `ImageDataDiff` objects have the following diff attributes:

    - ``diff_dimensions``: If the two arrays contain either a different number
      of dimensions or different sizes in any dimension, this contains a
      2-tuple of the shapes of each array.  Currently no further comparison is
      performed on images that don't have the exact same dimensions.

    - ``diff_pixels``: If the two images contain any different pixels, this
      contains a list of 2-tuples of the array index where the difference was
      found, and another 2-tuple containing the different values.  For example,
      if the pixel at (0, 0) contains different values this would look like::

          [(0, 0), (1.1, 2.2)]

      where 1.1 and 2.2 are the values of that pixel in each array.  This
      array only contains up to ``self.numdiffs`` differences, for storage
      efficiency.

    - ``diff_total``: The total number of different pixels found between the
      arrays.  Although ``diff_pixels`` does not necessarily contain all the
      different pixel values, this can be used to get a count of the total
      number of differences found.

    - ``diff_ratio``: Contains the ratio of ``diff_total`` to the total number
      of pixels in the arrays.
    """

    def __init__(self, a, b, numdiffs=10, rtol=0.0, atol=0.0, report_pixel_loc_diffs=False):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

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

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found.
        """

        self.report_pixel_loc_diffs = report_pixel_loc_diffs

        super().__init__(a, b,
                         numdiffs=numdiffs,
                         rtol=rtol,
                         atol=atol)

    def _diff(self):
        shapea, shapeb = self.a.shape, self.b.shape
        if shapea != shapeb:
            self.diff_dimensions = (shapea, shapeb)
            # Don't do any further comparison if the dimensions differ
            # TODO: Perhaps we could, however, diff just the intersection
            # between the two images
            return

        # If neither a nor b are floating point (or complex), ignore rtol and
        # atol
        if not (
                np.issubdtype(self.a.dtype, np.inexact)
                or np.issubdtype(self.b.dtype, np.inexact)
        ):
            rtol = 0
            atol = 0
        else:
            rtol = self.rtol
            atol = self.atol

        if self.report_pixel_loc_diffs:
            # Find the indices where the values are not equal
            diffs = where_not_allclose(self.a, self.b, atol=atol, rtol=rtol)

            self.diff_total = len(diffs[0])

            if self.diff_total == 0:
                # Then we're done
                return

            if self.numdiffs < 0:
                numdiffs = self.diff_total
            else:
                numdiffs = self.numdiffs

            self.diff_pixels = [
                (idx, (self.a[idx], self.b[idx]))
                for idx in islice(zip(*diffs), 0, numdiffs)
            ]
            self.diff_ratio = float(self.diff_total) / float(len(self.a.flat))

        else:
            # Make sure to separate nans in comparison
            data_within_tol = True

            nansa, nansb = np.isnan(self.a), np.isnan(self.b)
            a, b = self.a[~nansa], self.b[~nansb]
            # Only check the nans if the array values are the same
            if nansa.shape != nansb.shape or a.shape != b.shape:
                # Don't care about the actual numbers or locations, just set to something high
                data_within_tol = False
                self.diff_total = 999
                self.diff_ratio = 999.0
            elif a.shape == b.shape:
                # Check if data is within the tolerances (the non-nan data arrays are the same shape)
                if shapea == 4:
                    for nint in range(shapea[0]):
                        for ngrp in range(shapea[1]):
                            diff_total = np.abs(a[nint, ngrp, ...] - b[nint, ngrp, ...]) > (
                                    atol + rtol * np.abs(b[nint, ngrp, ...]))
                            self.diff_total = a[diff_total].size
                            if self.diff_total != 0:
                                data_within_tol = False
                                break
                        if not data_within_tol:
                            break
                if shapea == 3:
                    for ngrp in range(shapea[0]):
                        diff_total = np.abs(a[ngrp, ...] - b[ngrp, ...]) > (atol + rtol * np.abs(b[ngrp, ...]))
                        self.diff_total = a[diff_total].size
                        if self.diff_total != 0:
                            data_within_tol = False
                            break
                else:
                    diff_total = np.abs(a - b) > (atol + rtol * np.abs(b))
                    self.diff_total = a[diff_total].size
                    if self.diff_total != 0:
                        data_within_tol = False

            if not data_within_tol:
                # Don't care about the actual numbers or locations, just set to something high
                self.diff_ratio = 999.0
            else:
                # Data is the same, nothing to do
                self.diff_ratio = 0
                self.diff_total = 0
                return

    def _report(self):
        if self.diff_dimensions:
            dimsa = " x ".join(str(d) for d in reversed(self.diff_dimensions[0]))
            dimsb = " x ".join(str(d) for d in reversed(self.diff_dimensions[1]))
            self._writeln(" Data dimensions differ:")
            self._writeln(f"  a: {dimsa}")
            self._writeln(f"  b: {dimsb}")
            # For now, we don't do any further comparison if the dimensions
            # differ; though in the future it might be nice to be able to
            # compare at least where the images intersect
            self._writeln(" No further data comparison performed.")
            return

        if self.report_pixel_loc_diffs:

            if not self.diff_pixels:
                return

            max_relative = 0
            max_absolute = 0

            for index, values in self.diff_pixels:
                # Convert to int to avoid np.int64 in list repr.
                index = [int(x + 1) for x in reversed(index)]
                self._writeln(f" Data differs at {index}:")
                report_diff_values(
                    values[0],
                    values[1],
                    fileobj=self._fileobj,
                    indent_width=self._indent + 1,
                    rtol=self.rtol,
                    atol=self.atol,
                )
                rdiff, adiff = _get_differences(values[0], values[1])
                max_relative = max(max_relative, rdiff)
                max_absolute = max(max_absolute, adiff)

            if self.diff_total > self.numdiffs:
                self._writeln(" ...")
            self._writeln(
                f" {self.diff_total} different pixels found "
                f"({self.diff_ratio:.2%} different)."
            )
            self._writeln(f" Maximum relative difference: {max_relative}")
            self._writeln(f" Maximum absolute difference: {max_absolute}")


class STRawDataDiff(STImageDataDiff):
    """
    `RawDataDiff` is just a special case of `ImageDataDiff` where the images
    are one-dimensional, and the data is treated as a 1-dimensional array of
    bytes instead of pixel values.  This is used to compare the data of two
    non-standard extension HDUs that were not recognized as containing image or
    table data.

    `ImageDataDiff` objects have the following diff attributes:

    - ``diff_dimensions``: Same as the ``diff_dimensions`` attribute of
      `ImageDataDiff` objects. Though the "dimension" of each array is just an
      integer representing the number of bytes in the data.

    - ``diff_bytes``: Like the ``diff_pixels`` attribute of `ImageDataDiff`
      objects, but renamed to reflect the minor semantic difference that these
      are raw bytes and not pixel values.  Also the indices are integers
      instead of tuples.

    - ``diff_total`` and ``diff_ratio``: Same as `ImageDataDiff`.
    """

    def __init__(self, a, b, numdiffs=10, report_pixel_loc_diffs=False):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        numdiffs : int, optional
            The number of pixel/table values to output when reporting HDU data
            differences.  Though the count of differences is the same either
            way, this allows controlling the number of different values that
            are kept in memory or output.  If a negative value is given, then
            numdiffs is treated as unlimited (default: 10).

        report_pixel_loc_diffs : bool, optional
            As for ImageDiff, this will report all the locations where
            differences are found but instead of pixels is byte locations.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs

        super().__init__(a, b, numdiffs=numdiffs)

    def _diff(self):
        super()._diff()
        if self.diff_dimensions:
            self.diff_dimensions = (
                self.diff_dimensions[0][0],
                self.diff_dimensions[1][0],
            )

        self.diff_bytes = [(x[0], y) for x, y in self.diff_pixels]
        del self.diff_pixels

    def _report(self):
        if self.diff_dimensions:
            self._writeln(" Data sizes differ:")
            self._writeln(f"  a: {self.diff_dimensions[0]} bytes")
            self._writeln(f"  b: {self.diff_dimensions[1]} bytes")
            # For now we don't do any further comparison if the dimensions
            # differ; though in the future it might be nice to be able to
            # compare at least where the images intersect
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_bytes:
            return

        if self.report_pixel_loc_diffs:
            for index, values in self.diff_bytes:
                self._writeln(f" Data differs at byte {index}:")
                report_diff_values(
                    values[0],
                    values[1],
                    fileobj=self._fileobj,
                    indent_width=self._indent + 1,
                    rtol=self.rtol,
                    atol=self.atol,
                )

            self._writeln(" ...")
            self._writeln(
                f" {self.diff_total} different bytes found "
                f"({self.diff_ratio:.2%} different)."
            )


class STTableDataDiff(TableDataDiff):
    """
    Diff two table data arrays. It doesn't matter whether the data originally
    came from a binary or ASCII table--the data should be passed in as a
    recarray.

    `TableDataDiff` objects have the following diff attributes:

    - ``diff_column_count``: If the tables being compared have different
      numbers of columns, this contains a 2-tuple of the column count in each
      table.  Even if the tables have different column counts, an attempt is
      still made to compare any columns they have in common.

    - ``diff_columns``: If either table contains columns unique to that table,
      either in name or format, this contains a 2-tuple of lists. The first
      element is a list of columns (these are full `Column` objects) that
      appear only in table a.  The second element is a list of tables that
      appear only in table b.  This only lists columns with different column
      definitions, and has nothing to do with the data in those columns.

    - ``diff_column_names``: This is like ``diff_columns``, but lists only the
      names of columns unique to either table, rather than the full `Column`
      objects.

    - ``diff_column_attributes``: Lists columns that are in both tables but
      have different secondary attributes, such as TUNIT or TDISP.  The format
      is a list of 2-tuples: The first a tuple of the column name and the
      attribute, the second a tuple of the different values.

    - ``diff_values``: `TableDataDiff` compares the data in each table on a
      column-by-column basis.  If any different data is found, it is added to
      this list.  The format of this list is similar to the ``diff_pixels``
      attribute on `ImageDataDiff` objects, though the "index" consists of a
      (column_name, row) tuple.  For example::

          [('TARGET', 0), ('NGC1001', 'NGC1002')]

      shows that the tables contain different values in the 0-th row of the
      'TARGET' column.

    - ``diff_total`` and ``diff_ratio``: Same as `ImageDataDiff`.

    `TableDataDiff` objects also have a ``common_columns`` attribute that lists
    the `Column` objects for columns that are identical in both tables, and a
    ``common_column_names`` attribute which contains a set of the names of
    those columns.

    """

    def __init__(self, a, b, ignore_fields=[], numdiffs=10, rtol=0.0, atol=0.0,
                 report_pixel_loc_diffs=False):
        """
        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

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

        report_pixel_loc_diffs : bool, optional
            As for ImageDiff, this will report all the locations where
            differences are found but instead of pixels is column locations.

        """

        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.total_diff_per_col = {}

        super().__init__(a, b,
                         ignore_fields=ignore_fields,
                         numdiffs=numdiffs,
                         rtol=rtol,
                         atol=atol)

    def _diff(self):
        # Much of the code for comparing columns is similar to the code for
        # comparing headers--consider refactoring
        colsa = self.a.columns
        colsb = self.b.columns

        if len(colsa) != len(colsb):
            self.diff_column_count = (len(colsa), len(colsb))

        # Even if the number of columns are unequal, we still do comparison of
        # any common columns
        colsa = {c.name.lower(): c for c in colsa}
        colsb = {c.name.lower(): c for c in colsb}

        if "*" in self.ignore_fields:
            # If all columns are to be ignored, ignore any further differences
            # between the columns
            return

        # Keep the user's original ignore_fields list for reporting purposes,
        # but internally use a case-insensitive version
        ignore_fields = {f.lower() for f in self.ignore_fields}

        # It might be nice if there were a cleaner way to do this, but for now
        # it'll do
        for fieldname in ignore_fields:
            fieldname = fieldname.lower()
            if fieldname in colsa:
                del colsa[fieldname]
            if fieldname in colsb:
                del colsb[fieldname]

        colsa_set = set(colsa.values())
        colsb_set = set(colsb.values())
        self.common_columns = sorted(
            colsa_set.intersection(colsb_set), key=operator.attrgetter("name")
        )

        self.common_column_names = {col.name.lower() for col in self.common_columns}

        left_only_columns = {
            col.name.lower(): col for col in colsa_set.difference(colsb_set)
        }
        right_only_columns = {
            col.name.lower(): col for col in colsb_set.difference(colsa_set)
        }

        if left_only_columns or right_only_columns:
            self.diff_columns = (left_only_columns, right_only_columns)
            self.diff_column_names = ([], [])

        if left_only_columns:
            for col in self.a.columns:
                if col.name.lower() in left_only_columns:
                    self.diff_column_names[0].append(col.name)

        if right_only_columns:
            for col in self.b.columns:
                if col.name.lower() in right_only_columns:
                    self.diff_column_names[1].append(col.name)

        # If the tables have a different number of rows, we don't compare the
        # columns right now.
        # TODO: It might be nice to optionally compare the first n rows where n
        # is the minimum of the row counts between the two tables.
        if len(self.a) != len(self.b):
            self.diff_rows = (len(self.a), len(self.b))
            return

        # If the tables contain no rows there's no data to compare, so we're
        # done at this point. (See ticket #178)
        if len(self.a) == len(self.b) == 0:
            return

        # Like in the old fitsdiff, compare tables on a column by column basis
        # The difficulty here is that, while FITS column names are meant to be
        # case-insensitive, Astropy still allows, for the sake of flexibility,
        # two columns with the same name but different case.  When columns are
        # accessed in FITS tables, a case-sensitive is tried first, and failing
        # that a case-insensitive match is made.
        # It's conceivable that the same column could appear in both tables
        # being compared, but with different case.
        # Though it *may* lead to inconsistencies in these rare cases, this
        # just assumes that there are no duplicated column names in either
        # table, and that the column names can be treated case-insensitively.
        for col in self.common_columns:
            name_lower = col.name.lower()
            if name_lower in ignore_fields:
                continue

            cola = colsa[name_lower]
            colb = colsb[name_lower]

            for attr, _ in _COL_ATTRS:
                vala = getattr(cola, attr, None)
                valb = getattr(colb, attr, None)
                if diff_values(vala, valb):
                    self.diff_column_attributes.append(
                        ((col.name.upper(), attr), (vala, valb))
                    )

            arra = self.a[col.name]
            arrb = self.b[col.name]

            if np.issubdtype(arra.dtype, np.floating) and np.issubdtype(
                arrb.dtype, np.floating
            ):
                diffs = where_not_allclose(arra, arrb, rtol=self.rtol, atol=self.atol)
            elif "P" in col.format or "Q" in col.format:
                diffs = (
                    [
                        idx
                        for idx in range(len(arra))
                        if not np.allclose(
                            arra[idx], arrb[idx], rtol=self.rtol, atol=self.atol
                        )
                    ],
                )
            else:
                diffs = np.where(arra != arrb)

            self.diff_total += len(set(diffs[0]))

            # Find the total differences per column
            if not self.report_pixel_loc_diffs:
                if len(set(diffs[0])) > 0:
                    if col.name not in self.total_diff_per_col:
                        self.total_diff_per_col[col.name] = len(set(diffs[0]))
                    else:
                        self.total_diff_per_col[col.name] += len(set(diffs[0]))

            if self.numdiffs >= 0:
                if len(self.diff_values) >= self.numdiffs:
                    # Don't save any more diff values
                    continue

                # Add no more diff'd values than this
                max_diffs = self.numdiffs - len(self.diff_values)
            else:
                max_diffs = len(diffs[0])

            last_seen_idx = None
            for idx in islice(diffs[0], 0, max_diffs):
                if idx == last_seen_idx:
                    # Skip duplicate indices, which my occur when the column
                    # data contains multi-dimensional values; we're only
                    # interested in storing row-by-row differences
                    continue
                last_seen_idx = idx
                self.diff_values.append(((col.name, idx), (arra[idx], arrb[idx])))

        total_values = len(self.a) * len(self.a.dtype.fields)
        self.diff_ratio = float(self.diff_total) / float(total_values)

    def _report(self):
        if self.diff_column_count:
            self._writeln(" Tables have different number of columns:")
            self._writeln(f"  a: {self.diff_column_count[0]}")
            self._writeln(f"  b: {self.diff_column_count[1]}")

        if self.diff_column_names:
            # Show columns with names unique to either table
            for name in self.diff_column_names[0]:
                format = self.diff_columns[0][name.lower()].format
                self._writeln(f" Extra column {name} of format {format} in a")
            for name in self.diff_column_names[1]:
                format = self.diff_columns[1][name.lower()].format
                self._writeln(f" Extra column {name} of format {format} in b")

        col_attrs = dict(_COL_ATTRS)
        # Now go through each table again and show columns with common
        # names but other property differences...
        for col_attr, vals in self.diff_column_attributes:
            name, attr = col_attr
            self._writeln(f" Column {name} has different {col_attrs[attr]}:")
            report_diff_values(
                vals[0],
                vals[1],
                fileobj=self._fileobj,
                indent_width=self._indent + 1,
                rtol=self.rtol,
                atol=self.atol,
            )

        if self.diff_rows:
            self._writeln(" Table rows differ:")
            self._writeln(f"  a: {self.diff_rows[0]}")
            self._writeln(f"  b: {self.diff_rows[1]}")
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_values:
            return

        if self.report_pixel_loc_diffs:
            # Finally, let's go through and report column data differences:
            for indx, values in self.diff_values:
                self._writeln(" Column {} data differs in row {}:".format(*indx))
                report_diff_values(
                    values[0],
                    values[1],
                    fileobj=self._fileobj,
                    indent_width=self._indent + 1,
                    rtol=self.rtol,
                    atol=self.atol,
                )

            if self.diff_values and self.numdiffs < self.diff_total:
                self._writeln(
                    f" ...{self.diff_total - self.numdiffs} additional difference(s) found."
                )

            if self.diff_total > self.numdiffs:
                self._writeln(" ...")

        else:
            for colname in self.total_diff_per_col:
                self._writeln(" Column {} data differs on {} values".format(colname,
                                                                            self.total_diff_per_col[colname]))

        self._writeln(
            f" {self.diff_total} different table data element(s) found "
            f"({self.diff_ratio:.2%} different)."
        )


