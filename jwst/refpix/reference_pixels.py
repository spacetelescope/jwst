#  Module for handling Reference Pixels
#  Final CCWG Recommendation of 6/2013:
#
# The reference pixel correction for the NIR detectors should be done\
# immediately following the zero frame subtraction.  We recommend that
# the following steps be taken in order for each frame of each exposure,
# with the option to turn each one off *on a per-SCA basis*.  The
# parameters should eventually be optimized for each detector and instrument.
#
#   (1) For each amplifier, the sigma-clipped mean values for all odd and all
#       even columns of the horizontal (both top and bottom) reference pixels
#       should be calculated.  These values should then be subtracted from
#       every pixel in the corresponding odd/even columns.  There should be
#       an option to turn this step off, and replace with a single
#       sigma-clipped mean value for all horizontal reference pixels in
#       each amplifier.
#
#   (2) The vertical (both left and right) reference pixels should be smoothed
#       with an N-pixel wide boxcar convolution, where N may depend on detector
#       and instrument (adopt N=10 as a default).  The median value of the 8
#       smoothed reference pixels in each row should then be multiplied by a
#       gain factor of some value between 0 (which effectively turns off the
#       correction) and 1 (for full subtraction, should be the default), with
#       the exact value to be tuned for each detector and instrument.  Finally,
#       these smoothed and scaled values should be subtracted from every pixel
#       in the corresponding row.

#
#  Subarray processing added 7/2018
#
#  For NIR exposures, if the value of the meta.exposure.noutputs attribute is 1,
#  calculate the clipped means of odd and even columns
#  in detector coordinates.  Subtract the odd mean from the odd columns, and
#  the even mean from the even columns.  If there are no reference pixels in the
#  subarray, omit the refpix step.
#
#  If the value of meta.exposure.noutputs is 4, calculate odd and even reference
#  values for each amplifier separately, if available, and subtract those values
#  from their corresponding data sections.  Also use side reference pixels if
#  available.
#
#  For MIRI subarray exposures, omit the refpix step.

import numpy as np
from scipy import stats
import logging

from stdatamodels.jwst.datamodels import dqflags

from ..lib import reffile_utils

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#
# NIR Reference section dictionaries are zero indexed and specify the values
# to be used in the following slice:
# (rowstart: rowstop, colstart:colstop)
# The 'stop' values are one more than the actual final row or column, in
# accordance with how Python slices work

NIR_reference_sections = {'A': {'top': (2044, 2048, 0, 512),
                                'bottom': (0, 4, 0, 512),
                                'side': (0, 2048, 0, 4),
                                'data': (0, 2048, 0, 512)},
                          'B': {'top': (2044, 2048, 512, 1024),
                                'bottom': (0, 4, 512, 1024),
                                'data': (0, 2048, 512, 1024)},
                          'C': {'top': (2044, 2048, 1024, 1536),
                                'bottom': (0, 4, 1024, 1536),
                                'data': (0, 2048, 1024, 1536)},
                          'D': {'top': (2044, 2048, 1536, 2048),
                                'bottom': (0, 4, 1536, 2048),
                                'side': (0, 2048, 2044, 2048),
                                'data': (0, 2048, 1536, 2048)}
                          }

#
# MIR Reference section dictionaries are zero indexed and specify the values
# to be used in the following slice:
# name ('left' or 'right'): (rowstart, rowstop, column)
# except the 'data' entry:
# 'data': (rowstart, rowstop, colstart, colstop, stride)

MIR_reference_sections = {'A': {'left': (0, 1024, 0),
                                'right': (0, 1024, 1028),
                                'data': (0, 1024, 0, 1032, 4)},
                          'B': {'left': (0, 1024, 1),
                                'right': (0, 1024, 1029),
                                'data': (0, 1024, 1, 1032, 4)},
                          'C': {'left': (0, 1024, 2),
                                'right': (0, 1024, 1030),
                                'data': (0, 1024, 2, 1032, 4)},
                          'D': {'left': (0, 1024, 3),
                                'right': (0, 1024, 1031),
                                'data': (0, 1024, 3, 1032, 4)}
                          }

#
# Status returns

REFPIX_OK = 0
BAD_REFERENCE_PIXELS = 1
SUBARRAY_DOESNTFIT = 2
SUBARRAY_SKIPPED = 3


class Dataset():
    """Base Class to handle passing stuff from routine to routine

    Parameters:
    -----------

    input_model: data model object
        Science data model to be corrected

    is_subarray: boolean
        flag that shows whether the dataset was created from subarray
        data

    odd_even_columns: boolean
        flag that controls whether odd and even-numbered columns are
        processed separately (NIR only)

    use_side_ref_pixels: boolean
        flag the controls whether the side reference pixels are used in
        the correction (NIR only)

    side_smoothing_length: integer
        smoothing length the use in calculating the running median of
        the side reference pixels (NIR only)

    side_gain: float
        gain to use in applying the side reference pixel correction
        (NIR only)

    odd_even_rows: boolean
        flag that controls whether odd and even-numbered rows are handled
        separately (MIR only)

"""

    def __init__(self, input_model,
                 odd_even_columns,
                 use_side_ref_pixels,
                 side_smoothing_length,
                 side_gain,
                 odd_even_rows):

        if (input_model.meta.subarray.xstart is None or
                input_model.meta.subarray.ystart is None or
                input_model.meta.subarray.xsize is None or
                input_model.meta.subarray.ysize is None):
            raise ValueError('subarray metadata not found')
        self.input_model = input_model

        is_subarray = False
        if reffile_utils.is_subarray(input_model):
            is_subarray = True
        self.is_subarray = is_subarray

        self.zeroframe_proc = False

        (nints, ngroups, nrows, ncols) = input_model.data.shape
        self.nints = nints
        self.ngroups = ngroups
        self.nrows = nrows
        self.ncols = ncols
        self.full_shape = (nrows, ncols)
        self.detector = input_model.meta.instrument.detector
        self.noutputs = input_model.meta.exposure.noutputs

        self.xstart = input_model.meta.subarray.xstart
        self.ystart = input_model.meta.subarray.ystart
        self.xsize = input_model.meta.subarray.xsize
        self.ysize = input_model.meta.subarray.ysize

        self.colstart = self.xstart - 1
        self.colstop = self.colstart + self.xsize
        self.rowstart = self.ystart - 1
        self.rowstop = self.rowstart + self.ysize

        self.odd_even_columns = odd_even_columns
        self.use_side_ref_pixels = use_side_ref_pixels
        self.side_smoothing_length = side_smoothing_length
        self.side_gain = side_gain
        self.odd_even_rows = odd_even_rows
        self.bad_reference_pixels = False

        # Define temp array for processing every group
        self.pixeldq = self.get_pixeldq()
        self.group = None

    def sigma_clip(self, data, dq, low=3.0, high=3.0):
        """Wrap the scipy.stats.sigmaclip so that data with zero variance
        is handled cleanly

        Parameters:
        -----------

        data: NDArray
            Array of pixels to be sigma-clipped

        dq: NDArray
            DQ array for data

        low: float
            lower clipping boundary, in standard deviations from the mean (default=3.0)

        high: float
            upper clipping boundary, in standard deviations from the mean (default=3.0)

        Returns:
        --------

        mean: float
            clipped mean of data array

        """

        #
        # Only calculate the clipped mean for pixels that don't have the DO_NOT_USE
        # DQ bit set
        goodpixels = np.where(np.bitwise_and(dq, dqflags.pixel['DO_NOT_USE']) == 0)
        #
        # If there are no good pixels, return None
        if len(goodpixels[0]) == 0:
            return None
        #
        # scipy routine fails if the pixels all have exactly the same value
        if np.std(data[goodpixels], dtype=np.float64) != 0.0:
            clipped_ref, lowlim, uplim = stats.sigmaclip(data[goodpixels],
                                                         low, high)
            mean = clipped_ref.mean()
        else:
            mean = data[goodpixels].mean(dtype=np.float64)

        return mean

    def get_pixeldq(self):
        """Get the properly sized version of the pixeldq array from the
        input model.

        Parameters
        ----------
        None

        Returns
        -------
        pixeldq : NDArray
            numpy array for the pixeldq data with the full shape of the detector

        """
        if self.is_subarray:
            # deal with subarrays by embedding the pixeldq array in a full-sized
            # array with DO_NOT_USE and REFERENCE_PIXEL dqflags bit set where the
            # reference pixels live, except where the data are embedded
            if self.detector[:3] == 'MIR':
                fullrows = 1024
                fullcols = 1032
            else:
                fullrows = 2048
                fullcols = 2048
            self.full_shape = (fullrows, fullcols)
            pixeldq = np.zeros(self.full_shape, dtype=self.input_model.pixeldq.dtype)
            refpixdq_dontuse = dqflags.pixel['DO_NOT_USE'] | dqflags.pixel['REFERENCE_PIXEL']
            pixeldq[0:4, :] = refpixdq_dontuse
            pixeldq[fullrows - 4:fullrows, :] = refpixdq_dontuse
            pixeldq[4:fullrows - 4, 0:4] = refpixdq_dontuse
            pixeldq[4:fullrows - 4, fullcols - 4:fullcols] = refpixdq_dontuse
            pixeldq[self.rowstart:self.rowstop, self.colstart:self.colstop] = self.input_model.pixeldq.copy()
        else:
            pixeldq = self.input_model.pixeldq.copy()
        return pixeldq

    def get_group(self, integration, group):
        """Get a properly sized copy of the array for each group

        Parameters
        ----------
        integration : int
            Index of the integration from the input model from which to extract
            the group array

        group : int
            Index of the group, within the integration, from which to extract
            the group array
        """
        if self.group is None:
            self.group = np.zeros(self.full_shape, dtype=self.input_model.data.dtype)
        if self.is_subarray:
            self.group[self.rowstart:self.rowstop, self.colstart:self.colstop] = self.input_model.data[integration, group].copy()
        else:
            self.group[:, :] = self.input_model.data[integration, group].copy()

    def restore_group(self, integration, group):
        """Replace input model data with processed group array

        Parameters
        ----------
        integration : int
            Index of the integration from the input model which needs to be
            updated with the newly processed group array

        group : int
            Index of the group, within the integration, which needs to be
            updated with the newly processed group array

        """
        if self.is_subarray:
            self.input_model.data[integration, group] = self.group[self.rowstart:self.rowstop, self.colstart:self.colstop]
        else:
            self.input_model.data[integration, group] = self.group.copy()

    def log_parameters(self):
        """Print out the parameters that are valid for this type of data, and
        those that aren't

        Parameters
        ----------
        input_model : JWST datamodel
            Datamodel being processed

        Returns
        -------
        None

        """
        is_NIR = isinstance(self, NIRDataset)
        if is_NIR:
            if not self.is_subarray:
                log.info('NIR full frame data')
                log.info('The following parameters are valid for this mode:')
                log.info(f'use_side_ref_pixels = {self.use_side_ref_pixels}')
                log.info(f'odd_even_columns = {self.odd_even_columns}')
                log.info(f'side_smoothing_length = {self.side_smoothing_length}')
                log.info(f'side_gain = {self.side_gain}')
                log.info('The following parameter is not applicable and is ignored:')
                log.info(f'odd_even_rows = {self.odd_even_rows}')
            else:
                log.info('NIR subarray data')
                # Transform the pixeldq array from DMS to detector coords
                self.DMS_to_detector_dq()
                ngoodside = self.count_good_side_refpixels()
                ngoodtopbottom = self.count_good_top_bottom_refpixels()
                # Re-assign the pixeldq array since we transformed it to detector space
                # and we don't want to do it again
                self.pixeldq = self.get_pixeldq()
                is_4amp = False
                if self.noutputs == 4:
                    is_4amp = True
                if is_4amp:
                    log.info('4 readout amplifiers used')
                    if (ngoodside + ngoodtopbottom) == 0:
                        log.info('No valid reference pixels.  This step will have no effect')
                    else:
                        log.info('The following parameters are valid for this mode:')
                        if ngoodtopbottom > 0:
                            log.info(f'odd_even_columns = {self.odd_even_columns}')
                        if ngoodside > 0:
                            log.info(f'use_side_ref_pixels = {self.use_side_ref_pixels}')
                            log.info(f'side_smoothing_length = {self.side_smoothing_length}')
                            log.info(f'side_gain = {self.side_gain}')
                        log.info('The following parameters are not applicable and are ignored')
                        if ngoodtopbottom == 0:
                            log.info(f'odd_even_columns = {self.odd_even_columns}')
                        if ngoodside == 0:
                            log.info(f'use_side_ref_pixels = {self.use_side_ref_pixels}')
                            log.info(f'side_smoothing_length = {self.side_smoothing_length}')
                            log.info(f'side_gain = {self.side_gain}')
                        log.info(f'odd_even_rows = {self.odd_even_rows}')
                else:
                    log.info('Single readout amplifier used')
                    if ngoodtopbottom == 0:
                        log.info('No valid reference pixels.  This step wil have no effect')
                    else:
                        log.info('The following parameter is valid for this mode:')
                        log.info(f'odd_even_columns = {self.odd_even_columns}')
                        log.info('The following parameters are not applicable and are ignored:')
                        log.info(f'use_side_ref_pixels = {self.use_side_ref_pixels}')
                        log.info(f'side_smoothing_length = {self.side_smoothing_length}')
                        log.info(f'side_gain = {self.side_gain}')
                        log.info(f'odd_even_rows = {self.odd_even_rows}')
        else:
            if not self.is_subarray:
                log.info('MIRI full frame data')
                log.info('The following parameter is valid for this mode:')
                log.info(f'odd_even_rows = {self.odd_even_rows}')
                log.info('The following parameters are not applicable and are ignored:')
                log.info(f'use_side_ref_pixels = {self.use_side_ref_pixels}')
                log.info(f'odd_even_columns = {self.odd_even_columns}')
                log.info(f'side_smoothing_length = {self.side_smoothing_length}')
                log.info(f'side_gain = {self.side_gain}')
            else:
                log.info('MIRI subarray data')
                log.info('refpix processing skipped for this mode')

    def count_good_side_refpixels(self):
        donotuse = dqflags.pixel['DO_NOT_USE']
        ngood = 0
        for amplifier in 'AD':
            rowstart, rowstop, colstart, colstop = NIR_reference_sections[amplifier]['side']
            good = np.where(np.bitwise_and(self.pixeldq[rowstart:rowstop, colstart:colstop], donotuse) != donotuse)
            ngood += len(good[0])
        return ngood

    def count_good_top_bottom_refpixels(self):
        donotuse = dqflags.pixel['DO_NOT_USE']
        ngood = 0
        for edge in ['top', 'bottom']:
            for amplifier in 'ABCD':
                rowstart, rowstop, colstart, colstop = NIR_reference_sections[amplifier][edge]
                good = np.where(np.bitwise_and(self.pixeldq[rowstart:rowstop, colstart:colstop], donotuse) != donotuse)
                ngood += len(good[0])
        return ngood


class NIRDataset(Dataset):
    """Generic NIR detector Class.

    Parameters
    ----------

    input_model: data model object
        Science data model to be corrected

    is_subarray: boolean
        flag that shows whether the dataset was created from subarray
        data

    odd_even_columns: boolean
        flag that controls whether odd and even-numbered columns are
        processed separately

    use_side_ref_pixels: boolean
        flag the controls whether the side reference pixels are used in
        the correction

    side_smoothing_length: integer
        smoothing length the use in calculating the running median of
        the side reference pixels

    side_gain: float
        gain to use in applying the side reference pixel correction

    """

    def __init__(self, input_model,
                 odd_even_columns,
                 use_side_ref_pixels,
                 side_smoothing_length,
                 side_gain):

        super(NIRDataset, self).__init__(input_model,
                                         odd_even_columns,
                                         use_side_ref_pixels,
                                         side_smoothing_length,
                                         side_gain,
                                         odd_even_rows=False)

#
#  Even though the recommendation specifies calculating the mean of the
#  combined top and bottom reference sections, there's a good chance we
#  might want to calculate them separately
#
    def collect_odd_refpixels(self, group, amplifier, top_or_bottom):
        """Collect up the reference pixels corresponding to odd-numbered
        rows (first, third, fifth, etc, corresponding to even array indices)

        Parameters:
        -----------

        group: NDArray
            The group that is being processed

        amplifier: string ['A'|'B'|'C'|'D']
            String corresponding to the amplifier being processed

        top_or_bottom: string ['top'|'bottom']
            String corresponding to whether top or bottom reference pixels
            are bing processed

        Returns:

        oddref: NDArray
            Array containing all the odd reference pixels

        odddq: NDArray
            Array containing all the odd dq values for those reference pixels

        """
        rowstart, rowstop, colstart, colstop = \
            NIR_reference_sections[amplifier][top_or_bottom]

        oddref = group[rowstart:rowstop, colstart:colstop: 2]
        odddq = self.pixeldq[rowstart:rowstop, colstart:colstop: 2]
        return oddref, odddq

    def collect_even_refpixels(self, group, amplifier, top_or_bottom):
        """Collect up the reference pixels corresponding to even-numbered
        rows (second, fourth, sixth, etc, corresponding to odd array indices)

        Parameters:
        -----------

        group: NDArray
            The group that is being processed

        amplifier: string ['A'|'B'|'C'|'D']
            String corresponding to the amplifier being processed

        top_or_bottom: string ['top'|'bottom']
            String corresponding to whether top or bottom reference pixels
            are bing processed

        Returns:

        oddref: NDArray
            Array containing all the odd reference pixels

        odddq: NDArray
            Array containing all the odd dq values for those reference pixels

        """

        rowstart, rowstop, colstart, colstop = \
            NIR_reference_sections[amplifier][top_or_bottom]
        #
        # Even columns start on the second column
        colstart = colstart + 1
        evenref = group[rowstart:rowstop, colstart:colstop: 2]
        evendq = self.pixeldq[rowstart:rowstop, colstart:colstop: 2]
        return evenref, evendq

    def get_odd_refvalue(self, group, amplifier, top_or_bottom):
        """Calculate the clipped mean of the counts in the reference pixels
        in odd-numbered columns

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        top_or_bottom: string (['top'|'bottom'])
            Processing top or bottom reference pixels?

        Returns:
        --------

        odd: float
            Value of the clipped mean of the reference pixels in odd-numbered
            columns

        """

        ref, dq = self.collect_odd_refpixels(group, amplifier, top_or_bottom)
        odd = self.sigma_clip(ref, dq)
        return odd

    def get_even_refvalue(self, group, amplifier, top_or_bottom):
        """Calculate the clipped mean of the counts in the reference pixels
        in even-numbered columns

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        top_or_bottom: string (['top'|'bottom'])
            Processing top or bottom reference pixels?

        Returns:
        --------

        even: float
            Value of the clipped mean of the reference pixels in even-numbered
            columns

        """

        ref, dq = self.collect_even_refpixels(group, amplifier, top_or_bottom)
        even = self.sigma_clip(ref, dq)
        return even

    def get_amplifier_refvalue(self, group, amplifier, top_or_bottom):
        """Calculate the reference pixel mean for a given amplifier

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        top_or_bottom: string (['top'|'bottom'])
            Processing top or bottom reference pixels?

        Returns:
        --------

        Either:
            odd: float
                Value of the clipped mean of the reference pixels in odd-numbered
                columns

            even: float
                 Value of the clipped mean of the reference pixels in even-numbered
                 columns

        Or:
            mean: float
                Value of the clipped mean of the reference pixels in both odd-numbered
                and even-numbered columns
        """

        if self.odd_even_columns:
            odd = self.get_odd_refvalue(group, amplifier, top_or_bottom)
            even = self.get_even_refvalue(group, amplifier, top_or_bottom)
            if odd is None or even is None:
                self.bad_reference_pixels = True
            return odd, even
        else:
            rowstart, rowstop, colstart, colstop = \
                NIR_reference_sections[amplifier][top_or_bottom]
            ref = group[rowstart:rowstop, colstart:colstop]
            dq = self.pixeldq[rowstart:rowstop, colstart:colstop]
            mean = self.sigma_clip(ref, dq)
            if mean is None:
                self.bad_reference_pixels = True
            return mean

    def get_refvalues(self, group):
        """Get the reference pixel values for each amplifier, odd and even columns
        and top and bottom reference pixels

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        Returns:
        --------

        refpix: dictionary
            Dictionary containing the clipped mean of the reference pixels for
            each amplifier, odd and even columns (if selected, otherwise all columns)
            and top and bottom.

        """

        refpix = {}
        for amplifier in 'ABCD':
            refpix[amplifier] = {}
            refpix[amplifier]['odd'] = {}
            refpix[amplifier]['even'] = {}
            for top_bottom in ('top', 'bottom'):
                refvalues = self.get_amplifier_refvalue(group, amplifier,
                                                        top_bottom)
                if self.odd_even_columns:
                    refpix[amplifier]['odd'][top_bottom] = refvalues[0]
                    refpix[amplifier]['even'][top_bottom] = refvalues[1]
                else:
                    refpix[amplifier][top_bottom] = refvalues
        return refpix

    def do_top_bottom_correction(self, group, refvalues):
        """Do the top/bottom correction

        Parameters:
        ----------

        group: NDArray
            Group that is being processed

        refvalues: dictionary
            Dictionary of reference pixel clipped means

        Returns:
        --------

        None

        Side Effect:
        ------------

        The parameter _group_ is corrected for the bias drift using the
        top and bottom reference pixels

        """
        for amplifier in 'ABCD':
            datarowstart, datarowstop, datacolstart, datacolstop = \
                NIR_reference_sections[amplifier]['data']
            if self.odd_even_columns:
                oddreftop = refvalues[amplifier]['odd']['top']
                oddrefbottom = refvalues[amplifier]['odd']['bottom']
                evenreftop = refvalues[amplifier]['even']['top']
                evenrefbottom = refvalues[amplifier]['even']['bottom']
                #
                # For now, just average the top and bottom corrections
                oddrefsignal = self.average_with_None(oddreftop, oddrefbottom)
                evenrefsignal = self.average_with_None(evenreftop, evenrefbottom)
                if oddrefsignal is not None and evenrefsignal is not None:
                    oddslice = (slice(datarowstart, datarowstop, 1),
                                slice(datacolstart, datacolstop, 2))
                    evenslice = (slice(datarowstart, datarowstop, 1),
                                 slice(datacolstart + 1, datacolstop, 2))
                    group[oddslice] = group[oddslice] - oddrefsignal
                    group[evenslice] = group[evenslice] - evenrefsignal
                else:
                    pass
            else:
                reftop = refvalues[amplifier]['top']
                refbottom = refvalues[amplifier]['bottom']
                refsignal = self.average_with_None(reftop, refbottom)
                if refsignal is not None:
                    dataslice = (slice(datarowstart, datarowstop, 1),
                                 slice(datacolstart, datacolstop, 1))
                    group[dataslice] = group[dataslice] - refsignal
                else:
                    pass
        return

    def average_with_None(self, a, b):
        """Average two numbers.  If one is None, return the
        other.  If both are None, return None

        Parameters:
        -----------

        a, b:    Numbers or None

        Returns:
        --------

        result = Number or None
        """

        if a is None and b is None:
            return None

        if a is None:
            return b

        elif b is None:
            return a

        else:
            return 0.5 * (a + b)

    def create_reflected(self, data, smoothing_length):
        """Make an array bigger by extending it at the top and bottom by
        an amount equal to .5(smoothing length-1)
        (as the smoothing length will be odd)
        The extension is a reflection of the ends of the input array

        Parameters:
        -----------

        data: NDArray
            input data array

        smoothing_length: integer (should be odd, will be converted if not)
            smoothing length.  Amount by which the input array is extended is
            smoothing_length // 2 at the bottom and smoothing_length // 2 at
            the top

        Returns:
        --------

        reflected: NDArray
            array that has been extended at the top and bottom by reflecting the
            first and last few rows

        """

        nrows, ncols = data.shape
        if smoothing_length % 2 == 0:
            log.info("Smoothing length must be odd, adding 1")
            smoothing_length = smoothing_length + 1
        newheight = nrows + smoothing_length - 1
        reflected = np.zeros((newheight, ncols), dtype=data.dtype)
        bufsize = smoothing_length // 2
        reflected[bufsize:bufsize + nrows] = data[:]
        reflected[:bufsize] = data[bufsize:0:-1]
        reflected[-(bufsize):] = data[-2:-(bufsize + 2):-1]
        return reflected

    def median_filter(self, data, dq, smoothing_length):
        """Simple median filter.  Run a box of the same width as the data and
        height = smoothing_length.  Reflect the data at the top and bottom

        Parameters:
        -----------

        data: NDArray
            input 2-d science array

        dq: NDArray
            input 2-d dq array

        smoothing_length: integer (should be odd)
            height of box within which the median value is calculated

        Returns:
        --------

        result: NDArray
            1-d array that is a median filtered version of the input data
        """

        augmented_data = self.create_reflected(data, smoothing_length)
        augmented_dq = self.create_reflected(dq, smoothing_length)
        nrows, ncols = data.shape
        result = np.zeros(nrows)
        for i in range(nrows):
            rowstart = i
            rowstop = rowstart + smoothing_length
            goodpixels = np.where(np.bitwise_and(augmented_dq[rowstart:rowstop],
                                                 dqflags.pixel['DO_NOT_USE']) == 0)
            if len(goodpixels[0]) == 0:
                result[i] = np.nan
            else:
                window = augmented_data[rowstart:rowstop][goodpixels]
                result[i] = np.median(window)
        return result

    def calculate_side_ref_signal(self, group, colstart, colstop):
        """Calculate the reference pixel signal from the side reference pixels
        by running a box up the side reference pixels and calculating the running
        median

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        colstart: integer
            Starting column

        colstop: integer
            Ending column

        Returns:
        --------

        NDArray
            Median filtered version of the side reference pixels

        """

        smoothing_length = self.side_smoothing_length
        data = group[:, colstart:colstop + 1]
        dq = self.pixeldq[:, colstart:colstop + 1]
        return self.median_filter(data, dq, smoothing_length)

    def combine_ref_signals(self, left, right):
        """Combine the left and right reference signals by averaging
        on a row-by-row basis

        Parameters:
        -----------

        left: NDArray
            1-d array of median-filtered reference pixel values from the left side

        right: NDArray
            1-d array of median-filtered reference pixel values from the right side

        Returns:
        --------

        sidegroup: NDArray
            2-d array of average reference pixel vector replicated horizontally

        """

        combined = self.combine_with_NaNs(left, right)
        sidegroup = np.zeros((2048, 2048))
        for column in range(2048):
            sidegroup[:, column] = combined
        return sidegroup

    def combine_with_NaNs(self, a, b):
        """Combine 2 1-d arrays that have NaNs.
        Wherever both arrays are NaN, output is 0.0.
        Wherever a is NaN and b is not, return b.
        Wherever b is NaN and a is not, return a.
        Wherever neither a nor b is NaN, return the average of
        a and b

        Parameters:
        -----------

        a, b:   numpy 1-d arrays of numbers

        Returns:

        result = numpy 1-d array of numbers
        """

        result = np.zeros(len(a), dtype=a.dtype)

        bothnan = np.where(np.isnan(a) & np.isnan(b))
        result[bothnan] = 0.0

        a_nan = np.where(np.isnan(a) & ~np.isnan(b))
        result[a_nan] = b[a_nan]

        b_nan = np.where(~np.isnan(a) & np.isnan(b))
        result[b_nan] = a[b_nan]

        no_nan = np.where(~np.isnan(a) & ~np.isnan(b))
        result[no_nan] = 0.5 * (a[no_nan] + b[no_nan])

        return result

    def apply_side_correction(self, group, sidegroup):
        """Apply reference pixel correction from the side reference pixels

        Parameters:
        -----------

        group: NDArray
            Group being processed

        sidegroup: NDArray
            Side reference pixel signal replicated horizontally

        Returns:
        --------

        corrected_group: NDArray
            The group corrected for the side reference pixel signal

        """

        corrected_group = group - self.side_gain * sidegroup
        return corrected_group

    def do_side_correction(self, group):
        """Do all the steps of the side reference pixel correction

        Parameters:
        -----------

        group: NDArray
            Group being processed

        Returns:
        --------

        corrected_group: NDArray
            Corrected group

        """

        left = self.calculate_side_ref_signal(group, 0, 3)
        right = self.calculate_side_ref_signal(group, 2044, 2047)
        sidegroup = self.combine_ref_signals(left, right)
        corrected_group = self.apply_side_correction(group, sidegroup)
        return corrected_group

    def do_corrections(self):
        if self.is_subarray:
            if self.noutputs == 4:
                self.do_fullframe_corrections()
            else:
                self.do_subarray_corrections()
        else:
            self.do_fullframe_corrections()

    def do_fullframe_corrections(self):
        """Do Reference Pixels Corrections for all amplifiers, NIR detectors
        First read of each integration is NOT subtracted, as the signal is removed
        in the superbias subtraction step"""
        #
        #  First transform pixeldq array to detector coordinates
        self.DMS_to_detector_dq()

        for integration in range(self.nints):
            for group in range(self.ngroups):
                #
                # Get the reference values from the top and bottom reference
                # pixels
                #
                self.DMS_to_detector(integration, group)
                thisgroup = self.group
                refvalues = self.get_refvalues(thisgroup)
                self.do_top_bottom_correction(thisgroup, refvalues)
                if self.use_side_ref_pixels:
                    corrected_group = self.do_side_correction(thisgroup)
                    self.group = corrected_group
                else:
                    self.group = thisgroup
                #
                #  Now transform back from detector to DMS coordinates.
                self.detector_to_DMS(integration, group)
        log.setLevel(logging.INFO)
        return

    def do_subarray_corrections(self):
        """Do corrections for subarray.  Reference pixel value calculated
        separately for odd and even columns if odd_even_columns is True,
        otherwise a single number calculated from all reference pixels"""
        #
        #  First transform to detector coordinates
        #
        refdq = dqflags.pixel['REFERENCE_PIXEL']
        donotuse = dqflags.pixel['DO_NOT_USE']
        #
        # This transforms the pixeldq array from DMS to detector coordinates,
        # only needs to be done once
        self.DMS_to_detector_dq()
        # Determined refpix indices to use on each group
        refpixindices = np.where((self.pixeldq & refdq == refdq) & (self.pixeldq & donotuse != donotuse))
        nrefpixels = len(refpixindices[0])
        if nrefpixels == 0:
            self.bad_reference_pixels = True
            return
        if self.odd_even_columns:
            oddrefpixindices_row = []
            oddrefpixindices_col = []
            evenrefpixindices_row = []
            evenrefpixindices_col = []
            for i in range(nrefpixels):
                if (refpixindices[1][i] % 2) == 0:
                    evenrefpixindices_row.append(refpixindices[0][i])
                    evenrefpixindices_col.append(refpixindices[1][i])
                else:
                    oddrefpixindices_row.append(refpixindices[0][i])
                    oddrefpixindices_col.append(refpixindices[1][i])
            evenrefpixindices = (np.array(evenrefpixindices_row),
                                 np.array(evenrefpixindices_col))
            oddrefpixindices = (np.array(oddrefpixindices_row),
                                np.array(oddrefpixindices_col))

        for integration in range(self.nints):
            for group in range(self.ngroups):
                #
                # Get the reference values from the top and bottom reference
                # pixels
                #
                self.DMS_to_detector(integration, group)
                thisgroup = self.group

                if self.odd_even_columns:
                    evenrefpixvalue = self.sigma_clip(thisgroup[evenrefpixindices],
                                                      self.pixeldq[evenrefpixindices])
                    oddrefpixvalue = self.sigma_clip(thisgroup[oddrefpixindices],
                                                     self.pixeldq[oddrefpixindices])
                    thisgroup[:, 0::2] -= evenrefpixvalue
                    thisgroup[:, 1::2] -= oddrefpixvalue
                else:
                    refpixvalue = self.sigma_clip(thisgroup[refpixindices],
                                                  self.pixeldq[refpixindices])
                    thisgroup -= refpixvalue
                #
                #  Now transform back from detector to DMS coordinates.
                self.detector_to_DMS(integration, group)
        log.setLevel(logging.INFO)
        return


class NRS1Dataset(NIRDataset):
    """For NRS1 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRS1 is just flipped over the line X=Y
        self.get_group(integration, group)
        self.group = np.swapaxes(self.group, 0, 1)

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = np.swapaxes(self.group, 0, 1)
        self.restore_group(integration, group)

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = np.swapaxes(self.pixeldq, 0, 1)


class NRS2Dataset(NIRDataset):
    """NRS2 Data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRS2 is flipped over the line Y=X, then rotated 180 degrees
        self.get_group(integration, group)
        self.group = np.swapaxes(self.group, 0, 1)[::-1, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = np.swapaxes(self.pixeldq, 0, 1)[::-1, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # The inverse is to rotate 180 degrees, then flip over the line Y=X
        self.group = np.swapaxes(self.group[::-1, ::-1], 0, 1)
        self.restore_group(integration, group)


class NRCA1Dataset(NIRDataset):
    """For NRCA1 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCA1 is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)


class NRCA2Dataset(NIRDataset):
    """For NRCA2 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCA2 is just flipped in Y
        self.get_group(integration, group)
        self.group = self.group[::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1]
        self.restore_group(integration, group)


class NRCA3Dataset(NIRDataset):
    """For NRCA3 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCA3 is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)


class NRCA4Dataset(NIRDataset):
    """For NRCA4 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCA4 is just flipped in Y
        self.get_group(integration, group)
        self.group = self.group[::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1]
        self.restore_group(integration, group)


class NRCALONGDataset(NIRDataset):
    """For NRCALONG data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCALONG is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)


class NRCB1Dataset(NIRDataset):
    """For NRCB1 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCB1 is just flipped in Y
        self.get_group(integration, group)
        self.group = self.group[::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1]
        self.restore_group(integration, group)


class NRCB2Dataset(NIRDataset):
    """For NRCB2 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCB2 is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)
        # self.pixeldq = self.pixeldq[:, ::-1]


class NRCB3Dataset(NIRDataset):
    """For NRCB3 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCB3 is just flipped in Y
        self.get_group(integration, group)
        self.group = self.group[::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1]
        self.restore_group(integration, group)


class NRCB4Dataset(NIRDataset):
    """For NRCB4 data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCB4 is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)


class NRCBLONGDataset(NIRDataset):
    """For NRCBLONG data"""

    def DMS_to_detector(self, integration, group):
        #
        # NRCBLONG is just flipped in Y
        self.get_group(integration, group)
        self.group = self.group[::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1]
        self.restore_group(integration, group)


class NIRISSDataset(NIRDataset):
    """For NIRISS data"""

    def DMS_to_detector(self, integration, group):
        #
        # NIRISS has a 180 degree rotation followed by a flip across the line
        # X=Y
        self.get_group(integration, group)
        self.group = np.swapaxes(self.group[::-1, ::-1], 0, 1)

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = np.swapaxes(self.pixeldq[::-1, ::-1], 0, 1)

    def detector_to_DMS(self, integration, group):
        #
        # Just flip and rotate back
        self.group = np.swapaxes(self.group, 0, 1)[::-1, ::-1]
        self.restore_group(integration, group)


class GUIDER1Dataset(NIRDataset):
    """For GUIDER1 data"""

    def DMS_to_detector(self, integration, group):
        #
        # GUIDER1 is flipped in X and Y
        self.get_group(integration, group)
        self.group = self.group[::-1, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[::-1, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[::-1, ::-1]
        self.restore_group(integration, group)


class GUIDER2Dataset(NIRDataset):
    """For GUIDER2 data"""

    def DMS_to_detector(self, integration, group):
        #
        # GUIDER2 is just flipped in X
        self.get_group(integration, group)
        self.group = self.group[:, ::-1]

    def DMS_to_detector_dq(self):
        # pixeldq only has to be done once
        self.pixeldq = self.pixeldq[:, ::-1]

    def detector_to_DMS(self, integration, group):
        #
        # Just flip back
        self.group = self.group[:, ::-1]
        self.restore_group(integration, group)


class MIRIDataset(Dataset):
    """For MIRI data

    Parameters:
    -----------

    input_model: data model object
        Science data model to be corrected

    is_subarray: boolean
        flag that shows whether the dataset was created from subarray
        data

    odd_even_rows: boolean
        Flag that controls whether odd and even-numbered rows are
        handled separately

    """

    def __init__(self, input_model,
                 odd_even_rows):

        super(MIRIDataset, self).__init__(input_model,
                                          odd_even_columns=False,
                                          use_side_ref_pixels=False,
                                          side_smoothing_length=False,
                                          side_gain=False,
                                          odd_even_rows=odd_even_rows)

    def DMS_to_detector(self, integration, group):
        #
        # MIRI data doesn't need transforming
        pass

    def detector_to_DMS(self, integration, group):
        #
        # Do the opposite of above
        pass

    def collect_odd_refpixels(self, group, amplifier, left_or_right):
        """Collect reference pixels from odd-numbered rows

        Parameters:
        -----------

        group: NDArray
            Group being processed

        amplifier: string
            Amplifier being processed (['A'|'B'|'C'|'D'])

        left_or_right: string
            Process left or right side reference pixels (['left'|'right'])

        Returns:
        --------

        oddref: NDArray
            Reference pixels from odd-numbered rows

        odddq: NDArray
            DQ values for reference pixels from odd-numbered rows

        """

        rowstart, rowstop, column = MIR_reference_sections[amplifier][left_or_right]
        oddref = group[rowstart:rowstop:2, column]
        odddq = self.pixeldq[rowstart:rowstop:2, column]
        return oddref, odddq

    def collect_even_refpixels(self, group, amplifier, left_or_right):
        """Collect reference pixels from even-numbered rows

        Parameters:
        -----------

        group: NDArray
            Group being processed

        amplifier: string
            Amplifier being processed (['A'|'B'|'C'|'D'])

        left_or_right: string
            Process left or right side reference pixels (['left'|'right'])

        Returns:
        --------

        evenref: NDArray
            Reference pixels from even-numbered rows

        evendq: NDArray
            DQ values for reference pixels from even-numbered rows

        """

        rowstart, rowstop, column = MIR_reference_sections[amplifier][left_or_right]
        #
        # Even reference pixels start on the second row
        rowstart = rowstart + 1
        evenref = group[rowstart:rowstop:2, column]
        evendq = self.pixeldq[rowstart:rowstop:2, column]
        return evenref, evendq

    def get_odd_refvalue(self, group, amplifier, left_or_right):
        """Calculate the clipped mean of the counts in the reference pixels
        in odd-numbered rows

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        left_or_right: string (['left'|'right'])
            Processing left or right reference pixels?

        Returns:
        --------

        odd: float
            Value of the clipped mean of the reference pixels in odd-numbered
            rows

        """

        ref, dq = self.collect_odd_refpixels(group, amplifier, left_or_right)
        odd = self.sigma_clip(ref, dq)
        return odd

    def get_even_refvalue(self, group, amplifier, left_or_right):
        """Calculate the clipped mean of the counts in the reference pixels
        in even-numbered rows

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        left_or_right: string (['left'|'right'])
            Processing left or right reference pixels?

        Returns:
        --------

        even: float
            Value of the clipped mean of the reference pixels in even-numbered
            rows

        """

        ref, dq = self.collect_even_refpixels(group, amplifier, left_or_right)
        even = self.sigma_clip(ref, dq)
        return even

    def get_amplifier_refvalue(self, group, amplifier, left_or_right):
        """Calculate the reference pixel mean for a given amplifier

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        amplifier: string (['A'|'B'|'C'|'D'])
            Amplifier that is being processed

        left_or_right: string (['left'|'right'])
            Processing left or right side reference pixels?

        Returns:
        --------

        Either:
            odd: float
                Value of the clipped mean of the reference pixels in odd-numbered
                rows

            even: float
                 Value of the clipped mean of the reference pixels in even-numbered
                 rows

        Or:
            mean: float
                Value of the clipped mean of the reference pixels in both odd-numbered
                and even-numbered rows
        """

        if self.odd_even_rows:
            odd = self.get_odd_refvalue(group, amplifier, left_or_right)
            even = self.get_even_refvalue(group, amplifier, left_or_right)
            if odd is None:
                log.warning("Odd rows for amplifier {} have no good reference pixels".format(amplifier))
                self.bad_reference_piels = True
            elif even is None:
                log.warning("Even rows for amplifier {} have no good reference pixels".format(amplifier))
                self.bad_reference_pixels = True
            return odd, even
        else:
            rowstart, rowstop, column = MIR_reference_sections[amplifier][left_or_right]
            ref = group[rowstart:rowstop, column]
            dq = self.pixeldq[rowstart:rowstop, column]
            mean = self.sigma_clip(ref, dq)
            if mean is None:
                self.bad_reference_pixels = True
            return mean

    def get_refvalues(self, group):
        """Get the reference pixel values for each amplifier, odd and even rows
        and left and right side reference pixels

        Parameters:
        -----------

        group: NDArray
            Group that is being processed

        Returns:
        --------

        refpix: dictionary
            Dictionary containing the clipped mean of the reference pixels for
            each amplifier, odd and even rows (if selected, otherwise all rows)
            and left and right.

        """

        refpix = {}
        for amplifier in 'ABCD':
            refpix[amplifier] = {}
            refpix[amplifier]['odd'] = {}
            refpix[amplifier]['even'] = {}
            for left_right in ('left', 'right'):
                refvalues = self.get_amplifier_refvalue(group, amplifier,
                                                        left_right)
                if self.odd_even_rows:
                    refpix[amplifier]['odd'][left_right] = refvalues[0]
                    refpix[amplifier]['even'][left_right] = refvalues[1]
                else:
                    refpix[amplifier][left_right] = refvalues
        return refpix

    def do_left_right_correction(self, group, refvalues):
        """Do the reference pixel correction

        Parameters:
        ----------

        group: NDArray
            Group that is being processed

        refvalues: dictionary
            Dictionary of reference pixel clipped means

        Returns:
        --------

        None

        Side Effect:
        ------------

        The parameter _group_ is corrected for the bias drift using the
        left and right side reference pixels

        """

        for amplifier in 'ABCD':
            datarowstart, datarowstop, datacolstart, datacolstop, stride = \
                MIR_reference_sections[amplifier]['data']
            if self.odd_even_rows:
                oddrefleft = refvalues[amplifier]['odd']['left']
                oddrefright = refvalues[amplifier]['odd']['right']
                evenrefleft = refvalues[amplifier]['even']['left']
                evenrefright = refvalues[amplifier]['even']['right']
                #
                # For now, just average the left and right corrections
                oddrefsignal = 0.5 * (oddrefleft + oddrefright)
                evenrefsignal = 0.5 * (evenrefleft + evenrefright)
                oddslice = (slice(datarowstart, datarowstop, 2),
                            slice(datacolstart, datacolstop, 4))
                evenslice = (slice(datarowstart + 1, datarowstop, 2),
                             slice(datacolstart, datacolstop, 4))
                group[oddslice] = group[oddslice] - oddrefsignal
                group[evenslice] = group[evenslice] - evenrefsignal
            else:
                refleft = refvalues[amplifier]['left']
                refright = refvalues[amplifier]['right']
                refsignal = 0.5 * (refleft + refright)
                dataslice = (slice(datarowstart, datarowstop, 1),
                             slice(datacolstart, datacolstop, 4))
                group[dataslice] = group[dataslice] - refsignal
        return

    def do_corrections(self):
        if self.is_subarray:
            self.do_subarray_corrections()
        else:
            self.do_fullframe_corrections()

    def do_subarray_corrections(self):
        log.warning("Refpix correction skipped for MIRI subarray")
        return

    def do_fullframe_corrections(self):
        """Do Reference Pixels Corrections for all amplifiers, MIRI detectors"""
        #
        #  First we need to subtract the first read of each integration

        first_read = np.zeros((self.nints, self.nrows, self.ncols))
        log.info('Subtracting initial read from each integration')

        for i in range(self.nints):
            first_read[i] = self.input_model.data[i, 0].copy()
            self.input_model.data[i] = self.input_model.data[i] - first_read[i]

        #
        #  First transform to detector coordinates
        #
        for integration in range(self.nints):
            #
            #  Don't process the first group as it's all zeros and the clipped
            #  mean will return NaN
            #
            for group in range(1, self.ngroups):
                #
                # Get the reference values from the top and bottom reference
                # pixels
                #
                self.get_group(integration, group)
                thisgroup = self.group
                refvalues = self.get_refvalues(thisgroup)
                if self.bad_reference_pixels:
                    log.warning("Group {} has no reference pixels".format(group))
                    break
                self.do_left_right_correction(thisgroup, refvalues)
                #
                #  Now transform back from detector to DMS coordinates and transfer results to output
                self.restore_group(integration, group)
        log.setLevel(logging.INFO)
        #
        #  All done, now add the first read back in
        log.info('Adding initial read back in')

        for i in range(self.nints):
            self.input_model.data[i] += first_read[i]

        del first_read
        return


def create_dataset(input_model,
                   odd_even_columns,
                   use_side_ref_pixels,
                   side_smoothing_length,
                   side_gain,
                   odd_even_rows):
    """Create a dataset object from an input model.

    Parameters:
    -----------

    input_model: data model object
        Science data model to be corrected

    odd_even_columns: boolean
        flag that controls whether odd and even-numbered columns are
        processed separately (NIR only)

    use_side_ref_pixels: boolean
        flag the controls whether the side reference pixels are used in
        the correction (NIR only)

    side_smoothing_length: integer
        smoothing length the use in calculating the running median of
        the side reference pixels (NIR only)

    side_gain: float
        gain to use in applying the side reference pixel correction
        (NIR only)

    odd_even_rows: boolean
        flag that controls whether odd and even-numbered rows are handled
        separately (MIR only)

    """
    detector = input_model.meta.instrument.detector

    if reffile_utils.is_subarray(input_model):
        colstart = input_model.meta.subarray.xstart - 1
        colstop = colstart + input_model.meta.subarray.xsize
        rowstart = input_model.meta.subarray.ystart - 1
        rowstop = rowstart + input_model.meta.subarray.ysize
        if rowstart < 0 or colstart < 0 \
           or rowstop > 2048 or colstop > 2048:
            return None

    if detector[:3] == 'MIR':
        return MIRIDataset(input_model,
                           odd_even_rows)
    elif detector == 'NRS1':
        return NRS1Dataset(input_model,
                           odd_even_columns,
                           use_side_ref_pixels,
                           side_smoothing_length,
                           side_gain)
    elif detector == 'NRS2':
        return NRS2Dataset(input_model,
                           odd_even_columns,
                           use_side_ref_pixels,
                           side_smoothing_length,
                           side_gain)
    elif detector == 'NRCA1':
        return NRCA1Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCA2':
        return NRCA2Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCA3':
        return NRCA3Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCA4':
        return NRCA4Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCALONG':
        return NRCALONGDataset(input_model,
                               odd_even_columns,
                               use_side_ref_pixels,
                               side_smoothing_length,
                               side_gain)
    elif detector == 'NRCB1':
        return NRCB1Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCB2':
        return NRCB2Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCB3':
        return NRCB3Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCB4':
        return NRCB4Dataset(input_model,
                            odd_even_columns,
                            use_side_ref_pixels,
                            side_smoothing_length,
                            side_gain)
    elif detector == 'NRCBLONG':
        return NRCBLONGDataset(input_model,
                               odd_even_columns,
                               use_side_ref_pixels,
                               side_smoothing_length,
                               side_gain)
    elif detector == 'NIS':
        return NIRISSDataset(input_model,
                             odd_even_columns,
                             use_side_ref_pixels,
                             side_smoothing_length,
                             side_gain)
    elif detector == 'GUIDER1':
        return GUIDER1Dataset(input_model,
                              odd_even_columns,
                              use_side_ref_pixels,
                              side_smoothing_length,
                              side_gain)
    elif detector == 'GUIDER2':
        return GUIDER2Dataset(input_model,
                              odd_even_columns,
                              use_side_ref_pixels,
                              side_smoothing_length,
                              side_gain)
    else:
        log.error('Unrecognized detector')
        return NIRDataset(input_model,
                          odd_even_columns,
                          use_side_ref_pixels,
                          side_smoothing_length,
                          side_gain)


def correct_model(input_model, odd_even_columns,
                  use_side_ref_pixels,
                  side_smoothing_length, side_gain,
                  odd_even_rows):
    """Wrapper to do Reference Pixel Correction on a JWST Model.
    Performs the correction on the datamodel

    Parameters:
    -----------

    input_model: jwst.datamodels.model
        Model to be corrected

    odd_even_columns: boolean
        flag that controls whether odd and even-numbered columns are
        processed separately (NIR only)

    use_side_ref_pixels: boolean
        flag the controls whether the side reference pixels are used in
        the correction (NIR only)

    side_smoothing_length: integer
        smoothing length the use in calculating the running median of
        the side reference pixels (NIR only)

    side_gain: float
        gain to use in applying the side reference pixel correction
        (NIR only)

    odd_even_rows: boolean
        flag that controls whether odd and even-numbered rows are handled
        separately (MIR only)

    """
    if input_model.meta.instrument.name == 'MIRI':
        if reffile_utils.is_subarray(input_model):
            log.warning("Refpix correction skipped for MIRI subarrays")
            return SUBARRAY_SKIPPED

    input_dataset = create_dataset(input_model,
                                   odd_even_columns,
                                   use_side_ref_pixels,
                                   side_smoothing_length,
                                   side_gain,
                                   odd_even_rows)

    if input_dataset is None:
        status = SUBARRAY_DOESNTFIT
        return status
    input_dataset.log_parameters()
    reference_pixel_correction(input_dataset)

    return REFPIX_OK


def reference_pixel_correction(input_dataset):
    """
    Do the Reference Pixel Correction.

    Parameters:
    -----------

    input_dataset: Dataset
        Dataset to be corrected

    Returns:
    --------

    input_dataset: Dataset
        Corrected dataset

    """

    input_dataset.do_corrections()

    if input_dataset.input_model.meta.exposure.zero_frame:
        process_zeroframe_correction(input_dataset)

    return


def process_zeroframe_correction(input_dataset):
    """
    Do the Reference Pixel Correction for the ZEROFRAME array.

    Parameters
    ----------

    input_dataset : Dataset
        Dataset to be corrected

    Returns
    -------

    input_dataset : Dataset
        Corrected dataset
    """
    # Setup input model for ZEROFRAME
    saved_values = save_science_values(input_dataset)
    setup_dataset_for_zeroframe(input_dataset, saved_values)

    # Run refpix correction on ZEROFRAME
    input_dataset.do_corrections()

    restore_input_model(input_dataset, saved_values)


def restore_input_model(input_dataset, saved_values):
    """
    Restore the input model with saved values and move
    the computed ZEROFRAME value to the correct class
    variable.

    Parameters
    ----------
    input_dataset : Dataset
        Dataset to be corrected

    saved_values : tuple
        A tuple of saved values to be used to setup the final
        corrected RampModel.
    """
    data, gdq, pdq, wh_zero = saved_values

    nints, ngroups, nrows, ncols = data.shape
    zdims = (nints, nrows, ncols)

    # Get ZEROFRAME data
    zframe = input_dataset.input_model.data
    zdq = input_dataset.input_model.groupdq

    # Restore SCI data
    input_dataset.input_model.data = data
    input_dataset.input_model.groupdq = gdq
    input_dataset.input_model.pixeldq = pdq

    # Save computed ZEROFRAME
    zframe[zdq != 0] = 0.
    input_dataset.input_model.zeroframe = zframe.reshape(zdims)
    input_dataset.input_model.zeroframe[wh_zero] = 0.


def setup_dataset_for_zeroframe(input_dataset, saved_values):
    """
    Saves off corrected data for the SCI data.

    Parameters:
    -----------
    input_dataset : Dataset
        Dataset to be corrected
    """
    # Setup dimensions
    dims = input_dataset.input_model.zeroframe.shape
    nints, nrows, ncols = dims
    ngroups = 1
    new_dims = (nints, ngroups, nrows, ncols)

    # Setup ZEROFRAME data
    data = input_dataset.input_model.zeroframe
    data = data.reshape(new_dims)

    # Setup ZEROFRAME dummy groupdq
    gdtype = input_dataset.input_model.groupdq.dtype
    gdq = np.zeros(dims, dtype=gdtype)
    wh_zero = saved_values[-1]
    gdq[wh_zero] = dqflags.pixel['DO_NOT_USE']
    gdq = gdq.reshape(new_dims)

    # Setup dataset with ZEROFRAME data
    input_dataset.ngroups = ngroups
    input_dataset.pixeldq = input_dataset.get_pixeldq()
    input_dataset.input_model.data = data
    input_dataset.input_model.groupdq = gdq
    input_dataset.zeroframe_proc = True


def save_science_values(input_dataset):
    """
    Saves off corrected data for the SCI data.

    Parameters:
    -----------
    input_dataset : Dataset
        Dataset to be corrected

    Returns
    -------
    data : ndarray
        The correct SCI data.

    gdq : ndarray
        The correct SCI groupdq.

    pdq : ndarray
        The correct SCI pixeldq.

    wh_zero : ndarray
        The location of the zeroed out locations in the ZEROFRAME.
    """
    data = input_dataset.input_model.data
    gdq = input_dataset.input_model.groupdq
    pdq = input_dataset.input_model.pixeldq
    wh_zero = np.where(input_dataset.input_model.zeroframe[:, :, :] == 0.)

    return data, gdq, pdq, wh_zero
