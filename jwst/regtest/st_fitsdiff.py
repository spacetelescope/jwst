# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""STScI edits to astropy FitsDiff."""

import fnmatch
import operator
import textwrap
import warnings
from itertools import islice

import numpy as np
from astropy import __version__ as astropy_version
from astropy.io.fits.diff import (
    _COL_ATTRS,
    FITSDiff,
    HDUDiff,
    HeaderDiff,
    ImageDataDiff,
    TableDataDiff,
)
from astropy.io.fits.hdu.table import _TableLikeHDU
from astropy.table import Table
from astropy.utils.diff import diff_values, report_diff_values, where_not_allclose
from packaging.version import Version

__all__ = [
    "STFITSDiff",
    "STHDUDiff",
    "HeaderDiff",
    "STImageDataDiff",
    "STRawDataDiff",
    "STTableDataDiff",
]

#
# When we get to the right version, the following items need to be removed:
# Items                                                         Script where items live
# ----------------------------------------------------------    -----------------------
# global variable ASTROPY_LT_7_1_1                              st_fitsdiff.py
# function where_not_allclose                                   st_fitsdiff.py
# lines in function report_to_list that pass the abs/rel max    test_stfitsdiff.py
ASTROPY_LT_7_1_1 = Version(astropy_version) < Version("7.1.1.dev")

if ASTROPY_LT_7_1_1:
    # This function is for now copied from astropy. We can remove it once jwst updates the
    # astropy pin that contains this new functionality, i.e. above v.7.1.0.
    def where_not_allclose(a, b, rtol=1e-5, atol=1e-8, return_maxdiff=False):  # noqa: no-redef
        """
        Return array where values are above tolerances. Include max and min.

        A version of :func:`numpy.allclose` that returns the indices
        where the two arrays differ, instead of just a boolean value.

        Parameters
        ----------
        a, b : array-like
            Input arrays to compare.
        rtol, atol : float
            Relative and absolute tolerances as accepted by
            :func:`numpy.allclose`.
        return_maxdiff : bool
            Return the maximum of absolute and relative differences.

        Returns
        -------
        idx : tuple of array
            Indices where the two arrays differ.
        max_absolute : float
            Maximum of absolute difference, returned if ``return_maxdiff=True``.
        max_relative : float
            Maximum of relative difference, returned if ``return_maxdiff=True``.
        """
        # Create fixed mask arrays to handle INF and NaN; currently INF and NaN
        # are handled as equivalent
        a = np.ma.masked_invalid(a)
        b = np.ma.masked_invalid(b)

        absolute = np.ma.abs(b - a)

        if atol == 0.0 and rtol == 0.0:
            # Use a faster comparison for the most simple (and common) case
            thresh = 0
        else:
            thresh = atol + rtol * np.abs(b)

        # values invalid in only one of the two arrays should be reported
        invalid = a.mask ^ b.mask
        indices = np.where(invalid | (absolute.filled(0) > thresh))

        if return_maxdiff:
            absolute[invalid] = np.ma.masked
            finites = ~absolute.mask
            absolute = absolute.compressed()
            if len(indices[0]) == 0 or absolute.size == 0:
                max_absolute = max_relative = 0
            else:
                # remove all invalid values before computing max differences
                relative = absolute / np.abs(b[finites])
                max_absolute = float(np.max(absolute))
                max_relative = np.max(relative)
            return indices, max_absolute, max_relative
        else:
            return indices


def set_variable_to_empty_list(variable):
    if variable is None:
        variable = []
    return variable


class STFITSDiff(FITSDiff):
    """FITSDiff class that just filters warnings from astropy FITSDiff."""

    def _report(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            super()._report()


class STFITSDiffBeta(FITSDiff):
    """
    FITSDiff class from astropy with STScI ad hoc changes for STScI regression test reports.

    STScI changes include making sure that the files to be compared contain the same extensions
    and versions, regardless of order; if not, it will be reported. An option was implemented
    to provide a relative and/or absolute tolerance per extension, either with name or number.
    Additionally, header-specific tolerances can be provided for the primary extension or an
    absolute and relative tolerance can be provided to use for all headers.

    Full documentation of the base class is provided at:
    https://docs.astropy.org/en/stable/io/fits/api/diff.html
    """

    def __init__(
        self,
        a,
        b,
        ignore_hdus=None,
        ignore_keywords=None,
        ignore_comments=None,
        ignore_fields=None,
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
        report_pixel_loc_diffs=False,
        extension_tolerances=None,
    ):
        """
        For full documentation on variables, see original astropy code.

        Parameters
        ----------
        a : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object.

        b : str or `HDUList`
            The filename of a FITS file on disk, or an `HDUList` object to
            compare to the first file.

        ignore_hdus : sequence, optional
            HDU names to ignore when comparing two FITS files or HDU lists.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers.

        ignore_comments : sequence, optional
            List of header keywords whose comments should be ignored.

        ignore_fields : sequence, optional
            Case-insensitive names of any table columns to ignore.

        numdiffs : int, optional
            Number of pixel/table values to output when reporting HDU data differences.

        rtol : float, optional
            Relative difference to allow when comparing two float values.

        atol : float, optional
            Allowed absolute difference when comparing two float values.

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored.

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain whitespace.

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found.

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
        self.original_hdu_order = []
        self.expected_extension_tolerances = {}
        if extension_tolerances is not None:
            # Make sure the given dict keys are all upper case
            for key in extension_tolerances:
                # Make sure both tolerances exist in key
                tols = extension_tolerances[key]
                if "rtol" not in tols:
                    tols["rtol"] = rtol
                if "atol" not in tols:
                    tols["atol"] = atol
                if isinstance(key, str):
                    self.expected_extension_tolerances[key.upper()] = tols
                else:
                    self.expected_extension_tolerances[key] = tols
            # Make sure the other extensions get a default relative and absolute tolerance
            if "DEFAULT" not in [str(key).upper() for key in extension_tolerances]:
                self.expected_extension_tolerances["DEFAULT"] = {"rtol": rtol, "atol": atol}
        # Check if there are any numbers in the HDUs to ignore, and put them in a separate list
        self.ignore_number_hdus = []
        if ignore_hdus is not None:
            for key in ignore_hdus:
                if isinstance(key, int):
                    self.ignore_number_hdus.append(key)
            # Only keep extension names, remove the numbers
            for key in self.ignore_number_hdus:
                ignore_hdus.remove(key)

        super().__init__(
            a,
            b,
            ignore_hdus=set_variable_to_empty_list(ignore_hdus),
            ignore_keywords=set_variable_to_empty_list(ignore_keywords),
            ignore_comments=set_variable_to_empty_list(ignore_comments),
            ignore_fields=set_variable_to_empty_list(ignore_fields),
            numdiffs=numdiffs,
            rtol=rtol,
            atol=atol,
            ignore_blanks=ignore_blanks,
            ignore_blank_cards=ignore_blank_cards,
        )

    def _diff(self):
        # The following lines are identical to the original FITSDiff code

        if len(self.a) != len(self.b):
            self.diff_hdu_count = (len(self.a), len(self.b))

        # Record filenames for use later in _report
        self.filenamea = self.a.filename()
        if not self.filenamea:
            self.filenamea = f"<{self.a.__class__.__name__} object at {id(self.a):#x}>"

        self.filenameb = self.b.filename()
        if not self.filenameb:
            self.filenameb = f"<{self.b.__class__.__name__} object at {id(self.b):#x}>"

        # The following lines are STScI's additions:
        # 1. Get a list of additional HDUs to ignore, when ignore_hdus contain wildcards
        #    (only in this case will self.ignore_hdu_patterns be populated).
        additional_hdus_to_ignore = []
        if self.ignore_hdu_patterns:
            a_names = [hdu.name for hdu in self.a]
            b_names = [hdu.name for hdu in self.b]
            for pattern in self.ignore_hdu_patterns:
                a_ignored = fnmatch.filter(a_names, pattern)
                b_ignored = fnmatch.filter(b_names, pattern)
                additional_hdus_to_ignore = set.intersection(set(a_ignored), set(b_ignored))

        # The following lines are STScI's additions:
        # 1. Make sure that the files to be compared contain the same extensions
        #    and versions, regardless of order; if not, it will be reported.
        # 2. Option to provide a relative and/or absolute tolerance per extension,
        #    either with name or number.

        # Check if the files contain the same extension names and find the intersection
        ext_namesa = [ext.name for ext in self.a]
        ext_namesb = [ext.name for ext in self.b]
        ext_not_in_both = set(ext_namesa).symmetric_difference(set(ext_namesb))
        ext_intersection = set.intersection(set(ext_namesa), set(ext_namesb))
        if self.diff_hdu_count:
            self.diff_extnames = (
                sorted(ext_namesa),
                sorted(ext_namesb),
                sorted(ext_intersection),
                sorted(ext_not_in_both),
            )

        # Set the tolerance for the headers
        if "HEADERS" in list(self.expected_extension_tolerances):
            self.header_tolerances["rtol"] = self.expected_extension_tolerances["HEADERS"]["rtol"]
            self.header_tolerances["atol"] = self.expected_extension_tolerances["HEADERS"]["atol"]

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
                    if self.expected_extension_tolerances:
                        if idxa in self.expected_extension_tolerances:
                            self.rtol = self.expected_extension_tolerances[idxa]["rtol"]
                            self.atol = self.expected_extension_tolerances[idxa]["atol"]
                        elif extname in self.expected_extension_tolerances:
                            self.rtol = self.expected_extension_tolerances[extname]["rtol"]
                            self.atol = self.expected_extension_tolerances[extname]["atol"]
                        else:
                            self.rtol = self.expected_extension_tolerances["DEFAULT"]["rtol"]
                            self.atol = self.expected_extension_tolerances["DEFAULT"]["atol"]

                    # Only do a comparison if the HDU is not in any list to be ignored
                    if (
                        extname not in self.ignore_hdus
                        and idxa not in self.ignore_number_hdus
                        and extname not in additional_hdus_to_ignore
                    ):
                        hdu_diff = STHDUDiff.fromdiff(self, self.a[idxa], self.b[idxb])
                        if not hdu_diff.identical:
                            self.diff_hdus.append((idxa, hdu_diff, extname, extver))
                            self.original_hdu_order.append(idxa)

    def _report(self):
        # The following lines are identical to the original FITSDiff code

        wrapper = textwrap.TextWrapper(initial_indent="  ", subsequent_indent="  ")

        self._fileobj.write("\n")
        self._writeln(f" fitsdiff: {astropy_version}")
        self._writeln(f" a: {self.filenamea}\n b: {self.filenameb}")

        if self.ignore_hdus:
            ignore_hdus = " ".join(sorted(self.ignore_hdus))
            self._writeln(" HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdus))

        if self.ignore_hdu_patterns:
            ignore_hdu_patterns = " ".join(sorted(self.ignore_hdu_patterns))
            self._writeln(" HDU(s) not to be compared:\n" + wrapper.fill(ignore_hdu_patterns))

        if self.ignore_keywords:
            ignore_keywords = " ".join(sorted(self.ignore_keywords))
            self._writeln(" Keyword(s) not to be compared:\n" + wrapper.fill(ignore_keywords))

        if self.ignore_comments:
            ignore_comments = " ".join(sorted(self.ignore_comments))
            self._writeln(
                " Keyword(s) whose comments are not to be compared:\n"
                + wrapper.fill(ignore_comments)
            )

        if self.ignore_fields:
            ignore_fields = " ".join(sorted(self.ignore_fields))
            self._writeln(" Table column(s) not to be compared:\n" + wrapper.fill(ignore_fields))

        self._writeln(f" Maximum number of different data values to be reported: {self.numdiffs}")

        # The following lines are STScI's additions:
        # 1. Report the tolerances used per extension
        # 2. Specify the extensions present in each file, in both, and missing

        if not self.expected_extension_tolerances:
            self._writeln(f"\n Relative tolerance: {self.rtol}, Absolute tolerance: {self.atol}")

        if self.diff_hdu_count:
            self._fileobj.write("\n")
            self._writeln("Files contain different numbers of HDUs:")
            self._writeln(f" a: {self.diff_hdu_count[0]}, {self.diff_extnames[0]}")
            self._writeln(f" b: {self.diff_hdu_count[1]}, {self.diff_extnames[1]}")
            self._writeln(f" Common HDUs: {list(self.diff_extnames[2])}")
            self._writeln(f" Missing HDUs: {list(self.diff_extnames[3])}")

            if not self.diff_hdus:
                self._fileobj.write("\n")
                self._writeln("No differences found between common HDUs.")
                return

        elif not self.diff_hdus:
            self._fileobj.write("\n")
            self._writeln("No differences found.")
            return

        self._fileobj.write("\n")
        for ix, orig_idx in enumerate(self.original_hdu_order):
            idx, hdu_diff, extname, extver = self.diff_hdus[ix]
            if not self.expected_extension_tolerances:
                if idx == 0:
                    self._fileobj.write("\n")
                    self._writeln("Primary HDU:")
                else:
                    self._fileobj.write("\n\n")
                    if extname:
                        self._writeln(f"Extension HDU {orig_idx} ({extname}, {extver}):")
                    else:
                        self._writeln(f"Extension HDU {orig_idx}:")
            else:
                self._fileobj.write("\n\n")
                self._writeln(f"Extension HDU {orig_idx} ({extname}, {extver}):")
                if ix in self.expected_extension_tolerances:
                    rtol = self.expected_extension_tolerances[ix]["rtol"]
                    atol = self.expected_extension_tolerances[ix]["atol"]
                elif extname in self.expected_extension_tolerances:
                    rtol = self.expected_extension_tolerances[extname]["rtol"]
                    atol = self.expected_extension_tolerances[extname]["atol"]
                else:
                    rtol = self.expected_extension_tolerances["DEFAULT"]["rtol"]
                    atol = self.expected_extension_tolerances["DEFAULT"]["atol"]
                self._writeln(f"\n  Relative tolerance: {rtol:.4g}, Absolute tolerance: {atol:.4g}")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                hdu_diff.report(self._fileobj, indent=self._indent + 1)


class STHDUDiff(HDUDiff):
    """
    HDUDiff class from astropy with the STScI ad hoc changes for STScI regression test reports.

    STScI changes include making sure to use the header- or extension-specific absolute and
    relative tolerances, and print a stats report, which includes the following information:
    - number of nans
    - number of no nan values
    - number of zeros
    - Percentage difference of a from b for 0.0 absolute tolerance and 0.1, 0.001, 0.0001,
      1e-4, 1e-5, 1e-6, 1e-7, and 0.0 for relative tolerance
    - mean value
    - maximum value
    - minimum value
    - mean absolute difference
    - standard deviation for absolute difference
    - mean relative difference
    - standard deviation for relative difference

    Full documentation of the class is provided at:
    https://docs.astropy.org/en/stable/io/fits/api/diff.html
    """

    def __init__(
        self,
        a,
        b,
        ignore_keywords=None,
        ignore_comments=None,
        ignore_fields=None,
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        ignore_blanks=True,
        ignore_blank_cards=True,
        report_pixel_loc_diffs=False,
        header_tolerances=None,
    ):
        """
        For full documentation on variables, see original astropy code.

        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        ignore_keywords : sequence, optional
            Header keywords to ignore when comparing two headers.

        ignore_comments : sequence, optional
            List of header keywords whose comments should be ignored.

        ignore_fields : sequence, optional
            Case-insensitive names of any table columns to ignore.

        numdiffs : int, optional
            Number of pixel/table values to output when reporting HDU data differences.

        rtol : float, optional
            Relative difference to allow when comparing two float values.

        atol : float, optional
            Allowed absolute difference when comparing two float values.

        ignore_blanks : bool, optional
            Ignore extra whitespace at the end of string values either in
            headers or data. Extra leading whitespace is not ignored.

        ignore_blank_cards : bool, optional
            Ignore all cards that are blank, i.e. they only contain whitespace.

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found.

        header_tolerances : dict, optional
            Dictionary with the relative and absolute tolerances for all headers.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        if header_tolerances is None:
            header_tolerances = {}
        self.header_tolerances = header_tolerances
        self.nans, self.percentages, self.stats = None, None, None
        self.diff_dimensions = ()

        super().__init__(
            a,
            b,
            ignore_keywords=set_variable_to_empty_list(ignore_keywords),
            ignore_comments=set_variable_to_empty_list(ignore_comments),
            ignore_fields=set_variable_to_empty_list(ignore_fields),
            numdiffs=numdiffs,
            rtol=rtol,
            atol=atol,
            ignore_blanks=ignore_blanks,
            ignore_blank_cards=ignore_blank_cards,
        )

    def _diff(self):
        # The following lines are identical to the original HDUDiff code

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

        # The following lines are STScI's additions:
        # 1. Make sure to use header-specific absolute and relative tolerances
        # 2. Generate the stats report, which includes the following
        #    - number of nans
        #    - number of no nan values
        #    - number of zeros
        #    - Percentage difference of a from b for 0.0 absolute tolerance and
        #      0.1, 0.001, 0.0001, 1e-4, 1e-5, 1e-6, 1e-7, and 0.0 for relative tolerance
        #    - mean value
        #    - maximum value
        #    - minimum value
        #    - mean absolute difference
        #    - standard deviation for absolute difference
        #    - mean relative difference
        #    - standard deviation for relative difference

        # If specific header tolerances were given set them temporarily
        if self.header_tolerances:
            rtol, atol = self.rtol, self.atol
            self.rtol, self.atol = self.header_tolerances["rtol"], self.header_tolerances["atol"]

        # Get the header differences
        self.diff_headers = HeaderDiff.fromdiff(self, self.a.header.copy(), self.b.header.copy())
        # Reset the object tolerances
        if self.header_tolerances:
            self.rtol, self.atol = rtol, atol

        def report_array_zero_size(arr, arrnans, arrzeros):
            if arrnans.size > 0:
                nans_in_arr = arrnans.size
                nonans_in_arr = arr.size - nans_in_arr
            else:
                nans_in_arr = 0
                nonans_in_arr = arr.size
            if arrzeros.size > 0:
                zeros_in_arr = arrzeros.size
            else:
                zeros_in_arr = 0
            return zeros_in_arr, nans_in_arr, nonans_in_arr

        def get_quick_report(a, b):
            report_zeros_nan = Table()
            report_zeros_nan["Quantity"] = [
                "zeros",
                "nans",
                "no-nans",
                "min",
                "max",
                "mean",
                "std_dev",
            ]
            # Catch the case when the images are all nans and report accordingly
            nansa, nansb = a[np.isnan(a)], b[np.isnan(b)]
            zerosa, zerosb = a[a == 0.0], b[b == 0.0]
            zeros_in_a, nans_in_a, nonans_in_a = report_array_zero_size(a, nansa, zerosa)
            zeros_in_b, nans_in_b, nonans_in_b = report_array_zero_size(b, nansb, zerosb)
            nonansa, nonansb = a[np.isfinite(a)], b[np.isfinite(b)]
            mina, maxa, meana, stdeva = "-", "-", "-", "-"
            if nonansa.size > 0:
                mina = f"{np.min(nonansa):.4g}"
                maxa = f"{np.max(nonansa):.4g}"
                meana = f"{np.mean(nonansa):.4g}"
                stdeva = f"{np.std(nonansa):.4g}"
            minb, maxb, meanb, stdevb = "-", "-", "-", "-"
            if nonansb.size > 0:
                minb = f"{np.min(nonansb):.4g}"
                maxb = f"{np.max(nonansb):.4g}"
                meanb = f"{np.mean(nonansb):.4g}"
                stdevb = f"{np.std(nonansb):.4g}"
            # Populate report table
            report_zeros_nan["a"] = [zeros_in_a, nans_in_a, nonans_in_a, mina, maxa, meana, stdeva]
            report_zeros_nan["b"] = [zeros_in_b, nans_in_b, nonans_in_b, minb, maxb, meanb, stdevb]
            # Match nans for all arrays and remove them for logical comparison
            percentages, stats = Table(), Table()
            # Get the number of NaN in each array and other info
            n_total = b.size
            finite_idx = np.isfinite(a) & np.isfinite(b)
            finite_diffs = np.abs(b[finite_idx] - a[finite_idx])
            # Nothing to report if all values are 0 and the number of nans is the same
            # This is a failsafe but this bit of code will likely never be used
            # because arrays were found to be identical
            if (finite_diffs == 0.0).all() and nans_in_a == nans_in_b:
                return None, None, None
            # Calculate stats for absolute and relative differences
            # Catch the all NaNs case
            stats["Quantity"] = ["max", "min", "mean", "std_dev"]
            if finite_diffs.size == 0:
                percentages["threshold"] = [0.0]
                percentages["abs_diff%"] = [100]
                percentages["rel_diff%"] = [100]
                stats["abs_diff"] = [np.nan, np.nan, np.nan, np.nan]
                return report_zeros_nan, percentages, stats
            stats["abs_diff"] = [
                np.max(finite_diffs),
                np.min(finite_diffs),
                np.mean(finite_diffs),
                np.std(finite_diffs),
            ]
            stats["abs_diff"].format = "1.4g"
            relative_values = finite_diffs[b[finite_idx] != 0.0] / np.abs(
                b[finite_idx & (b != 0.0)]
            )
            # Catch an empty sequence
            if relative_values.size == 0:
                stats["rel_diff"] = [np.nan, np.nan, np.nan, np.nan]
            else:
                stats["rel_diff"] = [
                    np.max(relative_values),
                    np.min(relative_values),
                    np.mean(relative_values),
                    np.std(relative_values),
                ]
                stats["rel_diff"].format = "1.4g"
            # Calculate difference percentages
            thresholds = [0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 0.0]
            percentages["threshold"] = thresholds
            percent_abs_list = []
            for threshold in thresholds:
                n = finite_diffs[finite_diffs > (self.atol + threshold)].size
                percent_abs = float(n / n_total) * 100
                percent_abs_list.append(f"{percent_abs:.4g}")
            percentages["abs_diff%"] = percent_abs_list
            if relative_values.size > 0:
                percent_rel_list = []
                for threshold in thresholds:
                    n = relative_values[relative_values > (threshold + self.rtol)].size
                    percent_rel = float(n / n_total) * 100
                    percent_rel_list.append(f"{percent_rel:.4g}")
                percentages["rel_diff%"] = percent_rel_list
            return report_zeros_nan, percentages, stats

        # Code below contains mixed original HDUDiff lines as well as STScI's
        # to include the new classes and the stats reporting changes.

        if self.a.data is None or self.b.data is None:
            # TODO: Perhaps have some means of marking this case
            pass
        elif self.a.is_image and self.b.is_image:
            self.diff_data = STImageDataDiff.fromdiff(self, self.a.data, self.b.data)
            if self.diff_data.diff_total > 0:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    self.nans, self.percentages, self.stats = get_quick_report(
                        self.a.data, self.b.data
                    )
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
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    self.nans, self.percentages, self.stats = get_quick_report(
                        self.a.data, self.b.data
                    )
            # Clean up references to (possibly) memmapped arrays, so they can
            # be closed by .close()
            self.diff_data.a = None
            self.diff_data.b = None

    def _report(self):
        # The following lines are identical to the original HDUDiff code

        if self.identical:
            self._writeln(" No differences found.")
        if self.diff_extension_types:
            self._writeln(
                " Extension types differ:\n  a: {}\n  b: {}".format(*self.diff_extension_types)
            )
        if self.diff_extnames:
            self._writeln(" Extension names differ:\n  a: {}\n  b: {}".format(*self.diff_extnames))
        if self.diff_extvers:
            self._writeln(
                " Extension versions differ:\n  a: {}\n  b: {}".format(*self.diff_extvers)
            )

        if self.diff_extlevels:
            self._writeln(
                " Extension levels differ:\n  a: {}\n  b: {}".format(*self.diff_extlevels)
            )

        if not self.diff_headers.identical:
            self._fileobj.write("\n")
            self._writeln(" Headers contain differences:")
            self.diff_headers.report(self._fileobj, indent=self._indent + 1)

        # The following lines include STScI's changes to report the stats
        # as calculated in the _diff function above

        if self.diff_dimensions:
            dimsa = " x ".join(str(d) for d in reversed(self.diff_dimensions[0]))
            dimsb = " x ".join(str(d) for d in reversed(self.diff_dimensions[1]))
            self._writeln(" Data dimensions differ:")
            self._writeln(f"  a: {dimsa}")
            self._writeln(f"  b: {dimsb}")

        def report_data_diff():
            # Show differences in zeros and nans between a and b
            self._writeln("Values in a and b")
            for tline in self.nans.pformat():
                self._writeln(tline)

            # Show the difference (a-b) stats
            self._writeln("\nDifference stats: abs(b - a) ")
            for tline in self.stats.pformat():
                self._writeln(tline)

            # Show percentage differences
            self._writeln("\nPercentages of difference above (tolerance + threshold) ")
            for tline in self.percentages.pformat():
                self._writeln(tline)

        if self.diff_data is not None and not self.diff_data.identical:
            self._fileobj.write("\n")
            if [self.nans, self.percentages, self.stats] != [None, None, None]:
                report_data_diff()
            self.diff_data.report(self._fileobj, indent=self._indent + 1)


class STImageDataDiff(ImageDataDiff):
    """
    ImageDataDiff class with STScI ad hoc changes for STScI regression test reports.

    STScI changes include storing the shapes in variables to use later, and if the
    variable report_pixel_loc_diffs is True, print report as the original ImageDiff
    astropy class, otherwise if report_pixel_loc_diffs is False, use the shape to
    search for differences per integration or group, and truncate the loop and improve
    performance in large data; generate the ad hoc stats report described in
    the _diff function of the STHDUDiff class.

    Full documentation of the class is provided at:
    https://docs.astropy.org/en/stable/io/fits/api/diff.html
    """

    def __init__(self, a, b, numdiffs=10, rtol=0.0, atol=0.0, report_pixel_loc_diffs=False):
        """
        For full documentation on variables, see original astropy code.

        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        numdiffs : int, optional
            Number of pixel/table values to output when reporting HDU data differences.

        rtol : float, optional
            Relative difference to allow when comparing two float values.

        atol : float, optional
            Allowed absolute difference when comparing two float values.

        report_pixel_loc_diffs : bool, optional
            Report all the pixel locations where differences are found.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs

        super().__init__(a, b, numdiffs=numdiffs, rtol=rtol, atol=atol)

    def _diff(self):
        # Code below contains mixed original ImageDiff lines as well as STScI's
        # to include stats reporting changes:
        # 1. Store the shapes in variables to use later
        # 2. If the report_pixel_loc_diffs is True, print report as original ImageDiff
        #    report described in the _diff function of STHDUDiff.
        #    If report_pixel_loc_diffs is False:
        #    - Use the shape to search for differences per integration or group, and
        #      truncate the loop and improve performance in large data.
        #    - Generate the ad hoc stats report described in the _diff function of
        #      the STHDUDiff class

        shapea, shapeb = self.a.shape, self.b.shape
        if shapea != shapeb:
            self.diff_dimensions = (shapea, shapeb)
            # Don't do any further comparison if the dimensions differ
            # TODO: Perhaps we could, however, diff just the intersection
            # between the two images
            return

        rtol = self.rtol
        atol = self.atol

        if self.report_pixel_loc_diffs:
            # If neither a nor b are floating point (or complex), ignore rtol and
            # atol
            if not (
                np.issubdtype(self.a.dtype, np.inexact) or np.issubdtype(self.b.dtype, np.inexact)
            ):
                rtol = 0
                atol = 0

            # Find the indices where the values are not equal
            diffs, self.max_absolute, self.max_relative = where_not_allclose(
                self.a, self.b, atol=atol, rtol=rtol, return_maxdiff=True
            )

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
                for idx in islice(zip(*diffs, strict=False), 0, numdiffs)
            ]
            self.diff_ratio = float(self.diff_total) / float(len(self.a.flat))

        else:
            # Make sure to separate nans in comparison
            data_within_tol = True

            # Catch the case when the images are all nans
            nansa, nansb = self.a[np.isnan(self.a)], self.b[np.isnan(self.b)]
            nonana, nonanb = self.a[np.isfinite(self.a)], self.b[np.isfinite(self.b)]
            # Only check the nans if the array shapes are the same
            if nansa.shape != nansb.shape or nonana.shape != nonanb.shape:
                # No need to continue, there are differences. Go to stats calculation.
                data_within_tol = False
            else:
                # Check if data is within the tolerances (the non-nan data are the same shape)
                # but make sure that the nans are removed at the same place for both arrays
                a, b = self.a, self.b
                if len(shapea) == 4:
                    for nint in range(shapea[0]):
                        for ngrp in range(shapea[1]):
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", RuntimeWarning)
                                # Code that might generate a RuntimeWarning
                                # this data set is weird, do nothing and report
                                diff_total = np.abs(b[nint, ngrp, ...] - a[nint, ngrp, ...]) > (
                                    atol + rtol * np.abs(b[nint, ngrp, ...])
                                )
                            if a[nint, ngrp][diff_total].size != 0:
                                data_within_tol = False
                                break
                        if not data_within_tol:
                            break
                elif len(shapea) == 3:
                    for ngrp in range(shapea[0]):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            # Code that might generate a RuntimeWarning
                            # this data set is weird, do nothing and report
                            diff_total = np.abs(b[ngrp, ...] - a[ngrp, ...]) > (
                                atol + rtol * np.abs(b[ngrp, ...])
                            )
                            pass
                        if a[ngrp][diff_total].size != 0:
                            data_within_tol = False
                            break
                else:
                    finite_idx = np.isfinite(self.a) & np.isfinite(self.b)
                    a, b = self.a[finite_idx], self.b[finite_idx]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        # Code that might generate a RuntimeWarning
                        # this data set is weird, do nothing and report
                        diff_total = np.abs(b - a) > (atol + rtol * np.abs(b))
                        pass

                    if a[diff_total].size != 0:
                        data_within_tol = False

            if not data_within_tol:
                # Don't care about the actual numbers or locations, just set to something high
                self.diff_ratio = 999.0
                self.diff_total = 999
            else:
                # Data is the same, nothing to do
                self.diff_ratio = 0
                self.diff_total = 0

    def _report(self):
        # Code below contains mixed original ImageDiff lines as well as STScI's
        # to include stats reporting changes:
        # 1. No change with respect to original code for different dimensions
        # 2. Only print original report if report_pixel_loc_diffs is True

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

            self._writeln("\n * Pixel indices below are 1-based.")
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

            if self.diff_total > self.numdiffs:
                self._writeln(" ...")
            self._writeln(
                f" {self.diff_total} different pixels found ({self.diff_ratio:.2%} different)."
            )
            self._writeln(f" Maximum relative difference: {self.max_relative}")
            self._writeln(f" Maximum absolute difference: {self.max_absolute}")


class STRawDataDiff(STImageDataDiff):
    """
    Special instance of the STImageDataDiff class.

    STScI changes are to only print original report from the astropy RawDataDiff
    class if report_pixel_loc_diffs is True. Otherwise, generate the ad hoc stats
    report described report in the _diff function of the STHDUDiff class.

    Full documentation of the RawDataDiff class is provided at:
    https://docs.astropy.org/en/stable/io/fits/api/diff.html
    """

    def __init__(self, a, b, numdiffs=10, report_pixel_loc_diffs=False):
        """
        For full documentation on variables, see original astropy code.

        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        numdiffs : int, optional
            Number of pixel/table values to output when reporting HDU data differences.

        report_pixel_loc_diffs : bool, optional
            As for ImageDiff, this will report all the locations where
            differences are found but instead of pixels is byte locations.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs

        super().__init__(a, b, numdiffs=numdiffs)

    def _diff(self):
        # This function is exactly the same as the original, but it needs
        # to remain here because RawDataDiff is a special case of ImageDataDiff,
        # and since this was changed to use our version, the code would not
        # have access to the original _diff function in the RawDataDiff class.

        super()._diff()
        if self.diff_dimensions:
            self.diff_dimensions = (
                self.diff_dimensions[0][0],
                self.diff_dimensions[1][0],
            )

        self.diff_bytes = [(x[0], y) for x, y in self.diff_pixels]
        del self.diff_pixels

    def _report(self):
        # Code below contains mixed original RawDataDiff lines as well as STScI's
        # to include stats reporting changes.
        # 1. No change with respect to original code for different dimensions
        # 2. Only print original report if report_pixel_loc_diffs is True

        if self.diff_dimensions:
            self._writeln(" Data sizes differ:")
            self._writeln(f"  a: {self.diff_dimensions[0]} bytes")
            self._writeln(f"  b: {self.diff_dimensions[1]} bytes")
            # For now, we don't do any further comparison if the dimensions
            # differ; though in the future it might be nice to be able to
            # compare at least where the images intersect
            self._writeln(" No further data comparison performed.")
            return

        if not self.diff_bytes:
            return

        if self.report_pixel_loc_diffs:
            self._writeln("\n * Pixel indices below are 1-based.")
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
                f" {self.diff_total} different bytes found ({self.diff_ratio:.4g}% different)."
            )


class STTableDataDiff(TableDataDiff):
    """
    TableDataDiff class with STScI ad hoc changes for STScI regression test reports.

    STScI changes are to only report with the original format of the TableDataDiff astropy
    class if report_pixel_loc_diffs is True. Otherwise, the only relevant issue is if there
    are any differences, find out how many there are regardless of where they come from
    in the column. The report will only be the total differences per column.

    Full documentation of the class is provided at:
    https://docs.astropy.org/en/stable/io/fits/api/diff.html
    """

    def __init__(
        self,
        a,
        b,
        ignore_fields=None,
        numdiffs=10,
        rtol=0.0,
        atol=0.0,
        report_pixel_loc_diffs=False,
    ):
        """
        For full documentation on variables, see original astropy code.

        Parameters
        ----------
        a : BaseHDU
            An HDU object.

        b : BaseHDU
            An HDU object to compare to the first HDU object.

        ignore_fields : sequence, optional
            Case-insensitive names of any table columns to ignore.

        numdiffs : int, optional
            Number of pixel/table values to output when reporting HDU data differences.

        rtol : float, optional
            Relative difference to allow when comparing two float values
            either in header values, image arrays, or table columns.

        atol : float, optional
            Allowed absolute difference.

        report_pixel_loc_diffs : bool, optional
            As for ImageDiff, this will report all the locations where
            differences are found but instead of pixels is column locations.
        """
        self.report_pixel_loc_diffs = report_pixel_loc_diffs
        self.rel_diffs = 0
        self.report_table = Table(
            names=(
                "col_name",
                "dtype",
                "abs_diffs",
                "abs_max",
                "abs_mean",
                "abs_std",
                "rel_diffs",
                "rel_max",
                "rel_mean",
                "rel_std",
            ),
            dtype=(
                "str",
                "str",
                "int32",
                "float64",
                "float64",
                "float64",
                "float64",
                "float64",
                "float64",
                "float64",
            ),
        )
        self.report_zeros_nan = Table(
            names=(
                "col_name",
                "zeros_a_b",
                "nan_a_b",
                "no-nan_a_b",
                "max_a_b",
                "min_a_b",
                "mean_a_b",
            ),
            dtype=(
                "str",
                "str",
                "str",
                "str",
                "str",
                "str",
                "str",
            ),
        )

        super().__init__(
            a,
            b,
            ignore_fields=set_variable_to_empty_list(ignore_fields),
            numdiffs=numdiffs,
            rtol=rtol,
            atol=atol,
        )

    def _diff(self):
        # The following lines are identical to the original TableDataDiff code

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

        left_only_columns = {col.name.lower(): col for col in colsa_set.difference(colsb_set)}
        right_only_columns = {col.name.lower(): col for col in colsb_set.difference(colsa_set)}

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
                    self.diff_column_attributes.append(((col.name.upper(), attr), (vala, valb)))

            arra = self.a[col.name]
            arrb = self.b[col.name]

            # The following lines include STScI's changes:
            # - Calculate the absolute and relative differences separately for ad hoc report
            # - If report_pixel_loc_diffs is False, just get the total differences,
            #   it is not important where they come from

            def get_stats_if_nonans(nonans):
                if nonans.size > 0:
                    nonansmax = np.max(nonans)
                    nonansmin = np.min(nonans)
                    nonansmean = np.mean(nonans)
                else:
                    nonansmax = np.nan
                    nonansmin = np.nan
                    nonansmean = np.nan
                return nonansmax, nonansmin, nonansmean

            # Calculate the absolute and relative differences
            get_stats = False
            reported = False
            if np.issubdtype(arra.dtype, np.floating) and np.issubdtype(arrb.dtype, np.floating):
                # Catch the case of a column full of nans
                nansa = arra[np.isnan(arra)]
                nansb = arrb[np.isnan(arrb)]
                if nansa.size != nansb.size:
                    nonansa = arra[np.isfinite(arra)]
                    nonansb = arrb[np.isfinite(arrb)]
                    arramax, arramin, arramean = get_stats_if_nonans(nonansa)
                    arrbmax, arrbmin, arrbmean = get_stats_if_nonans(nonansb)
                    finite_idx = np.isfinite(arra) & np.isfinite(arrb)
                    diffs = np.abs(arrb[finite_idx] - arra[finite_idx])
                    maxa, meana, stda, rdiff, maxr, meanr, stdr = 0, 0, 0, 0, 0, 0, 0
                    if diffs.size > 0:
                        maxa = np.max(diffs)
                        meana = np.mean(diffs)
                        stda = np.std(diffs)
                        finiteb = arrb[finite_idx]
                        rel_diffs = diffs[finiteb != 0.0] / finiteb[finiteb != 0.0]
                        rdiff = rel_diffs.size
                        if rdiff > 0:
                            maxr = np.max(rel_diffs)
                            meanr = np.mean(rel_diffs)
                            stdr = np.std(rel_diffs)

                    # Report the total number of zeros, nans, and no-nan values
                    self.report_zeros_nan.add_row(
                        (
                            col.name,
                            f"{arra[arra == 0.0].size} {arrb[arrb == 0.0].size}",
                            f"{nansa.size} {nansb.size}",
                            f"{nonansa.size} {nonansb.size}",
                            f"{arramax:>9.4g} {arrbmax:>9.4g}",
                            f"{arramin:>9.4g} {arrbmin:>9.4g}",
                            f"{arramean:>9.4g} {arrbmean:>9.4g}",
                        )
                    )
                    self.report_table.add_row(
                        (
                            col.name,
                            str(arra.dtype).replace(">", ""),
                            diffs.size,
                            maxa,
                            meana,
                            stda,
                            rdiff,
                            maxr,
                            meanr,
                            stdr,
                        )
                    )

                    if not self.report_pixel_loc_diffs:
                        self.diff_total += arrb.size
                        self.rel_diffs += np.nan
                    reported = True
                else:
                    nan_idx = np.isnan(arra) | np.isnan(arrb)
                    anonan = arra[~nan_idx]
                    bnonan = arrb[~nan_idx]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        diffs = np.abs(bnonan - anonan)
                    abs_diffs = diffs[diffs > self.atol].size
                    nozeros = (diffs != 0.0) & (bnonan != 0.0)
                    rel_values = diffs[nozeros] / np.abs(bnonan[nozeros])
                    rel_diffs = rel_values[rel_values > self.rtol].size
                    get_stats = True

            elif "P" in col.format or "Q" in col.format:
                diffs = (
                    [
                        idx
                        for idx in range(len(arra))
                        if not np.allclose(arra[idx], arrb[idx], rtol=self.rtol, atol=self.atol)
                    ],
                )
                reported = True
                self.report_zeros_nan.add_row(
                    (
                        col.name,
                        "- -",
                        "- -",
                        "- -",
                        "- -",
                        "- -",
                        "- -",
                    )
                )
                abs_diffs = np.array(diffs).size
                self.report_table.add_row(
                    (
                        col.name,
                        str(arra.dtype).replace(">", ""),
                        abs_diffs,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    )
                )

                if not self.report_pixel_loc_diffs:
                    self.diff_total += abs_diffs
                    self.rel_diffs += 0

            if not reported:
                if get_stats:
                    sum_diffs = np.any(diffs > (self.atol + self.rtol * np.abs(bnonan)))
                    if sum_diffs and len(rel_values) > 0:
                        # Report the total number of zeros, nans, and no-nan values
                        self.report_zeros_nan.add_row(
                            (
                                col.name,
                                f"{arra[arra == 0.0].size} {arrb[arrb == 0.0].size}",
                                f"{arra[np.isnan(arra)].size} {arrb[np.isnan(arrb)].size}",
                                f"{arra[~np.isnan(arra)].size} {arrb[~np.isnan(arrb)].size}",
                                f"{np.max(anonan):>9.4g} {np.max(bnonan):>9.4g}",
                                f"{np.min(anonan):>9.4g} {np.min(bnonan):>9.4g}",
                                f"{np.mean(anonan):>9.4g} {np.mean(bnonan):>9.4g}",
                            )
                        )

                        self.report_table.add_row(
                            (
                                col.name,
                                str(arra.dtype).replace(">", ""),
                                abs_diffs,
                                np.max(diffs),
                                np.mean(diffs),
                                np.std(diffs),
                                rel_diffs,
                                np.max(rel_values),
                                np.mean(rel_values),
                                np.std(rel_values),
                            )
                        )

                        if not self.report_pixel_loc_diffs:
                            self.diff_total += abs_diffs
                            self.rel_diffs += rel_diffs

                else:
                    diffs = arra[arra != arrb].size
                    if diffs > 0:
                        # Report the differences per column
                        self.report_table.add_row(
                            (
                                col.name,
                                str(arra.dtype).replace(">", ""),
                                diffs,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                            )
                        )
                        if not self.report_pixel_loc_diffs:
                            self.diff_total += diffs
                            self.rel_diffs = np.nan

            if self.report_pixel_loc_diffs:
                # The following lines are identical to the original TableDataDiff code

                if np.issubdtype(arra.dtype, np.floating) and np.issubdtype(
                    arrb.dtype, np.floating
                ):
                    diffs, self.max_absolite, self.max_relative = where_not_allclose(
                        arra, arrb, rtol=self.rtol, atol=self.atol, return_maxdiff=True
                    )
                elif "P" in col.format or "Q" in col.format:
                    diffs = (
                        [
                            idx
                            for idx in range(len(arra))
                            if not np.allclose(arra[idx], arrb[idx], rtol=self.rtol, atol=self.atol)
                        ],
                    )
                else:
                    diffs = np.where(arra != arrb)

                self.diff_total += len(set(diffs[0]))

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

        # Calculate the absolute difference
        total_values = len(self.a) * len(self.a.dtype.fields)
        # Calculate the absolute and relative difference percentages
        if self.report_pixel_loc_diffs:
            total_values = len(self.a) * len(self.a.dtype.fields)
            self.diff_ratio = float(self.diff_total) / float(total_values)
        self.stdiff_ratio = (float(self.diff_total) / float(total_values)) * 100
        self.diff_ratio_rel = (float(self.rel_diffs) / float(total_values)) * 100

    def _report(self):
        # The following lines are identical to the original TableDataDiff code

        if self.diff_column_count:
            self._writeln(" Tables have different number of columns:")
            self._writeln(f"  a: {self.diff_column_count[0]}")
            self._writeln(f"  b: {self.diff_column_count[1]}")

        if self.diff_column_names:
            # Show columns with names unique to either table
            for name in self.diff_column_names[0]:
                colformat = self.diff_columns[0][name.lower()].format
                self._writeln(f" Extra column {name} of format {colformat} in a")
            for name in self.diff_column_names[1]:
                colformat = self.diff_columns[1][name.lower()].format
                self._writeln(f" Extra column {name} of format {colformat} in b")

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

        if not self.diff_total:
            return

        # The following lines include STScI's changes:
        # - Only report the locations of the differences if report_pixel_loc_diffs is True,
        #   otherwise report the column and the total number of differences and the
        #   percentage absolute and relative differences.

        self._writeln(
            f"Found {self.diff_total} different table data element(s). "
            "Reporting percentages above respective tolerances: "
            f"\n    - absolute .... {self.stdiff_ratio:.4g}%"
        )
        if np.isnan(self.diff_ratio_rel):
            self._writeln(
                "\n * Unable to calculate relative differences and stats due to data types"
            )
        else:
            self._writeln(f"    - relative .... {self.diff_ratio_rel:.4g}%")

        # Print differences in zeros and nans per column
        self._writeln("\nValues in a and b")
        tlines = self.report_zeros_nan.pformat(max_width=-1)
        for tline in tlines:
            self._writeln(tline)

        # Print the difference (b - a) stats
        self._writeln("\nDifference stats: abs(b - a) ")
        # make sure the format is acceptable
        for colname in self.report_table.columns:
            if colname in ["col_name", "dtype"]:
                continue
            self.report_table[colname].format = ".4g"
        tlines = self.report_table.pformat(max_width=-1)
        for tline in tlines:
            self._writeln(tline)

        # Report of column differences from astropy
        if self.report_pixel_loc_diffs:
            # Finally, let's go through and report column data differences:
            self._writeln("\n * Pixel indices below are 1-based.")
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

            self._writeln(
                f" {self.diff_total} different table data element(s) found "
                f"({self.diff_ratio:.2%} different)."
            )
