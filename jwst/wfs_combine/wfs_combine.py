import logging

import numpy as np
from scipy.signal import convolve, correlate2d
from scipy.ndimage import gaussian_filter, center_of_mass
from scipy.interpolate import griddata

from .. import datamodels
from ..datamodels import dqflags


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]


class DataSet:
    """
    Two dithered input wavefront sensing images to be combined

    """

    def __init__(self, infile_1, infile_2, outfile, do_refine, flip_dithers, psf_size,
                 blur_size, n_size):
        """
        Short Summary
        -------------
        Assign models to input files.

        Parameters:
        -----------
        infile_1: string
            First input file
        infile_2: string
            Second input file
        outfile: string
            File for combined image
        do_refine: boolean
            True if refined offset calculation and application is to be made
        flip_dithers: boolean
            True will cause the dithers to align in pixel coordinates for different filters
        psf_size: float
            size of largest psf
        blur_size: float
            amount of smoothing to apply before finding the initial centroid
        n_size: int
            size of interpolation box
        """

        if outfile == "":
            log.error('No output product specified in the association table.')

        try:
            self.input_1 = datamodels.open(infile_1)
            self.input_2 = datamodels.open(infile_2)
        except IOError:
            log.error('Error creating a model from at least 1 of : %s %s',
                      infile_1, infile_2)

        self.do_refine = do_refine

        if self.input_1.data.shape != self.input_2.data.shape:
            log.error('Incompatible sizes for input files')

        log.info('Output file: %s', outfile)
        log.info('do_refine: %s', do_refine)
        log.info('flip_dithers: %s', flip_dithers)
        self.file1 = infile_1
        self.file2 = infile_2
        self.off_x = 0
        self.off_y = 0
        self.flt_off_x = 0
        self.flt_off_y = 0
        self.diff = np.ndarray(shape=(2 * psf_size + 1, 2 * psf_size + 1))
        self.flip_dithers = flip_dithers
        self.psf_size = psf_size
        self.blur_size = blur_size
        self.n_size = n_size

    def do_all(self):
        """
        Short Summary
        -------------
        Execute all tasks for Wave Front Sensing Combination

        Parameters
        ----------

        Returns
        -------
        new_model: data model object
            combined input file data
        """

        self.off_x, self.off_y = self.get_wcs_offsets()

        # If the shift in x is negative, switch the two images
        if self.off_x < 0 and self.flip_dithers:
            self.input_1.close()
            self.input_2.close()
            self.input_1 = datamodels.open(self.file2)
            self.input_2 = datamodels.open(self.file1)
            log.info('File 1 to combine: %s', self.file2)
            log.info('File 2 to combine: %s', self.file1)
        else:
            log.info('File 1 to combine: %s', self.file1)
            log.info('File 2 to combine: %s', self.file2)

        # Input SCI arrays may have nan's so replace with 0's to prevent
        # later annoyances (hopefully this can be removed later)
        self.input_1.data[np.isnan(self.input_1.data)] = 0.
        self.input_2.data[np.isnan(self.input_2.data)] = 0.
        im_1_a = self.input_1.copy()  # Aligned image #1 (already aligned)
        im_2_a = self.create_aligned_2()  # Aligned image #2

        # Create and populate extensions for combined data
        data_c, dq_c, err_c, diff_c = self.create_combined(im_1_a, im_2_a)

        log.info(f"Final x, y offset in pixels: {self.off_x} {self.off_y}")

        self.diff = diff_c
        # Create a new model using the combined arrays...
        new_model = datamodels.ImageModel(data=data_c, dq=dq_c, err=err_c)
        new_model.update(self.input_1)
        new_model.history.append('Flip dithers = {}'.format(self.flip_dithers))
        new_model.history.append('WFS_COMBINE refine offset = {}'.format(self.do_refine))
        new_model.history.append('WFS_COMBINE X offset applied ' + str(self.off_x) + ' pixels ' +
                                 'actual offset ' + str(round(self.flt_off_x, 2)) + ' pixels')
        new_model.history.append('WFS_COMBINE Y offset applied ' + str(self.off_y) + ' pixels ' +
                                 'actual offset ' + str(round(self.flt_off_y, 2)) + ' pixels')
        return new_model

    def create_aligned_2(self):
        """
        Short Summary
        -------------
        Align image 2 in image 1's frame.

        Long Summary
        ------------
        If refined offset determination is selected, do steps 1-7 else do 7 only
        1. Create a smoothed image of the input SCI data of image #1. First
           create an image to smooth by first setting SCI pixels with bad DQ
           values equal to the mean of the good pixels. Then smooth this
           'repaired' image using a gaussian kernel of size BLUR_SIZE.
        2. Find the approximate centroid of this PSF, by taking all the pixels
           in this smoothed image that exceed 50% of the maximum of the
           smoothed image, and taking the mean of the coordinates of these
           pixels.  The x- and y-mean defines the centroid location.
        3. Set the limits of the subarrays for interpolation by taking this
           centroid +/- psf_size and adding the BLUR_SIZE, taking the edges
           into account.
        4. Determine overlap of these interpolated images, and return nominally
           aligned, interpolated images.
        5. Around this nominal alignment, calculate refined (delta) offsets.
        6. Add the refined delta offsets to the nominal offsets.
        7. Do final alignment for original (not interpolated) image #2

        Parameters
        ----------

        Returns
        -------
        model_2_a:  Data Model object
           aligned model for input image #2
        """

        self.off_x, self.off_y = self.get_wcs_offsets()
        log.info(f"x,y offset in integer pixels from WCS: {self.off_x} {self.off_y}")
        if self.do_refine:
            # 1. Create smoothed image of input SCI data of image #1
            # 1a. create image to smooth by first setting bad DQ pixels equal
            #     to mean of good pixels
            data_1 = self.input_1.data.astype(float)
            bad1 = np.bitwise_and(self.input_1.dq, DO_NOT_USE).astype(bool)
            data_1[bad1] = data_1[~bad1].mean()

            # 1b. Create smoothed image by smoothing this 'repaired' image
            g = gauss_kern(self.blur_size, sizey=None)
            s_data_1 = convolve(data_1, g, mode='valid')

            # 2. Find approximate center of PSF in unsmoothed frame by taking
            #    all pixels in smoothed image exceeding 50% of the maximum
            #    of the smoothed image, and taking the mean of the coordinates
            #    of these pixels. Add BLUR_SIZE to take smoothing into account
            wh_data_hi = np.where(s_data_1 > 0.5 * s_data_1.max())

            ctrd_x = wh_data_hi[1].mean() + self.blur_size
            ctrd_y = wh_data_hi[0].mean() + self.blur_size

            log.info('Approximate centroid of image 1 PSF has x,y : %s %s',
                     round(ctrd_x), round(ctrd_y))

            # 3. Set limits of the subarrays (in frames of input data)
            #    for interpolation by taking this centroid +/- psf_size
            #    and adding BLUR_SIZE, taking edges into account
            xmin = int(round(max(0, ctrd_x - self.psf_size)))
            ymin = int(round(max(0, ctrd_y - self.psf_size)))
            xmax = int(round(min(self.input_1.data.shape[1], ctrd_x + self.psf_size)))
            ymax = int(round(min(self.input_1.data.shape[0], ctrd_y + self.psf_size)))

            # 3a. Set subarrays and interpolate over bad pixels
            data_sub_1 = self.input_1.data[ymin: ymax, xmin: xmax]
            dq_sub_1 = self.input_1.dq[ymin: ymax, xmin: xmax]
            sci_int_1 = interp_array(data_sub_1, dq_sub_1, self.n_size)

            data_sub_2 = self.input_2.data[ymin: ymax, xmin: xmax]
            dq_sub_2 = self.input_2.dq[ymin: ymax, xmin: xmax]
            sci_int_2 = interp_array(data_sub_2, dq_sub_2, self.n_size)

            # 4. Determine overlap of these interpolated images, and
            #    return nominally aligned, interpolated images
            sci_nai_1, sci_nai_2 = get_overlap(sci_int_1, sci_int_2,
                                               self.off_x, self.off_y)
            # 5. Around this nominal alignment, get refined (delta) offsets
            ref_del_off_x, ref_del_off_y = calc_refined_offsets(sci_nai_1, sci_nai_2, 0, 0, self.psf_size)
            log.info('From the refined offsets calculation,'
                     'the x,y changes in offsets are: %s %s',
                     round(ref_del_off_x, 2), round(ref_del_off_y, 2))

            # 6. Add the refined delta offsets to the nominal offsets
            self.flt_off_x = self.off_x + ref_del_off_x
            self.flt_off_y = self.off_y + ref_del_off_y
            self.off_x += int(round(ref_del_off_x))
            self.off_y += int(round(ref_del_off_y))

        # Do the final alignment for original (not interpolated) image two
        data_2_a, dq_2_a, err_2_a = self.apply_final_offsets()

        model_2_a = self.input_2.copy()  # Model for aligned image #2
        model_2_a.data = data_2_a
        model_2_a.dq = dq_2_a
        model_2_a.err = err_2_a

        return model_2_a

    def apply_final_offsets(self):
        """
        Short Summary
        -------------
        Apply final offsets, aligning each array for image #2 to #1's frame

        Parameters
        ----------

        Returns
        -------
        data_2_a: 2D float array
            aligned SCI array of image #2
        dq_2_a: 2D int array
            aligned DQ array of image #2
        err_2_a: 2D float array
            aligned ERR array of image #2
        """

        data_2_a = self.do_2d_shifts(self.input_2.data)
        dq_2_a = self.do_2d_shifts(self.input_2.dq)
        err_2_a = self.do_2d_shifts(self.input_2.err)

        return data_2_a, dq_2_a, err_2_a

    def get_wcs_offsets(self):
        """
        Short Summary
        -------------
        Get the nominal offsets from the WCS information of each of the
        2 input DataModel objects. From the difference in pointings (in
        pixels) of the 2 images, round off to the nearest integers as
        the specifications require that the pointings will differ by exact
        integers.

        Parameters
        ----------

        Returns
        -------
        off_x: integer
            difference (#2 -#1) in pointing in pixels in the x-direction
        off_y: integer
            difference (#2 -#1) in pointing in pixels in the y-direction
        """
        wcs1 = self.input_1.meta.wcs
        wcs2 = self.input_2.meta.wcs
        tr1 = wcs1.get_transform('detector', 'world')
        tr2 = wcs2.get_transform('world', 'detector')

        # Get coords of center pixel
        xcen = int(self.input_1.data.shape[1] / 2)
        ycen = int(self.input_1.data.shape[0] / 2)

        radec = tr1(xcen, ycen)
        pixels = tr2(radec[0], radec[1])
        off_x = pixels[0] - xcen
        off_y = pixels[1] - ycen
        log.info('From the WCS the x,y pixel offsets are: %s %s',
                 round(off_x, 2), round(off_y, 2))
        self.flt_off_x = off_x
        self.flt_off_y = off_y
        off_x = int(round(off_x))  # Offsets required to be integers
        off_y = int(round(off_y))

        return off_x, off_y

    def create_combined(self, image1, image2):
        """
        Short Summary
        -------------
        Create combined image from aligned input images. In the combined image:

        The SCI pixel values are set by:
        1. for pixels that are good (based on DQ) in both images, use their average
        2. for pixels that are good in image #1 and bad in image #2, use image #1
        3. for pixels that are bad in image #1 and good in image #2, use image #2
        4. for pixels that are bad in both images, leave as default (0)

        The DQ pixel values are set by:
        1. use pixels that are good in either image #1 or image #2
        2. for pixels that are bad in both images, add a 'DO_NOT_USE' value to the
           corresponding DQ value

        The ERR pixel values are similarly set:
        1. for pixels that are good in both images, use their average (will modify
           later)
        2. for pixels that are good in image #1 and bad in image #2, use image #1
        3. for pixels that are bad in image #1 and good in image #2, use image #2
        4. for pixels that are bad in both images, leave as default (0)

        The WCS of the output model is set to the WCS of the 1st input

        Parameters
        ----------
        image1: 2D Data Model
             aligned image from input #1
        image2: 2D Data Model
             aligned image from input #2

        Returns
        -------
        data_comb: 2d float array
            combined SCI array
        dq_comb: 2d integer array
            combined DQ array
        err_comb: 2d float array
            combined ERR array
        """

        data1 = image1.data.astype(float)
        data2 = image2.data.astype(float)
        dq1 = image1.dq.copy()
        dq2 = image2.dq.copy()
        err1 = image1.err.copy()
        err2 = image2.err.copy()

        # Create boolean arrays of bad pixels in each input image
        bad1 = np.bitwise_and(dq1, DO_NOT_USE).astype(bool)
        good1 = ~bad1
        bad2 = np.bitwise_and(dq2, DO_NOT_USE).astype(bool)
        good2 = ~bad2

        # Combine via algorithm set out above

        # Data pixels that are bad in both will stay 0
        data_comb = np.zeros_like(data1)
        data_comb[good1 & good2] = 0.5 * (data1[good1 & good2] + data2[good1 & good2])
        data_comb[good1 & bad2] = data1[good1 & bad2]
        data_comb[good2 & bad1] = data2[good2 & bad1]

        dq_comb = dq1.copy()
        dq_comb[good1 & bad2] = dq1[good1 & bad2]
        dq_comb[good2 & bad1] = dq2[good2 & bad1]
        dq_comb[bad1 & bad2] = np.bitwise_or(DO_NOT_USE, dq_comb[bad1 & bad2])

        err_comb = np.zeros_like(err1)
        err_comb[good1 & good2] = 0.5 * (err1[good1 & good2] + err2[good1 & good2])
        err_comb[good1 & bad2] = err1[good1 & bad2]
        err_comb[good2 & bad1] = err2[good2 & bad1]

        data_diff = np.zeros_like(data1)
        data_diff[good1 & good2] = (data1[good1 & good2] - data2[good1 & good2])
        data_diff[good1 & bad2] = 0
        data_diff[good2 & bad1] = 0

        return data_comb, dq_comb, err_comb, data_diff

    def do_2d_shifts(self, a):
        """
        Short Summary
        -------------
        Create 2d output array by shifting 2d array input by (off_x, off_y),
        where output has same dimensions as input.
        Parameters
        ----------
        a: 2d float array
            input array
        Returns
        -------
        b: 2d float array
            shifted array of input a
        """

        ai_x, af_x = get_final_index_range(self.off_x, a.shape[1])
        ai_y, af_y = get_final_index_range(self.off_y, a.shape[0])

        bi_x = a.shape[1] - af_x  # For output, x-direction's initial channel
        bf_x = a.shape[1] - ai_x  # ...and final channel

        bi_y = a.shape[0] - af_y  # For output, y-direction's initial channel
        bf_y = a.shape[0] - ai_y  # ...and final channel

        b = np.zeros(a.shape)
        b[bi_y:bf_y, bi_x:bf_x] = a[ai_y:af_y, ai_x:af_x]

        return b


def get_final_index_range(offset, length):
    """
    Short Summary
    -------------
    Get the initial and final indices for the given offset and array length:
    For offset <= 0: i = 0,  f = length - abs(offset)
    For offset > 0: i = offset,  f = length
    Parameters
    ----------
    offset: integer
        offset
    length: integer
        length of (1D) array
    Returns
    -------
    i: integer
        initial index
    f: integer
        final index
    """

    i = int((abs(offset) + offset) / 2)
    f = length + int((-abs(offset) + offset) / 2)

    return i, f


def gauss_kern(size, sizey=None):
    """
    Short Summary
    -------------
    Returns a normalized 2D gauss kernel array for convolution.

    Parameters
    ----------
    size: int
        size of gaussian kernel in x_dim
    sizey: int
        sizey of gaussian kernel in y_dim

    Returns
    -------
    g/g.sum(): 2D float array
        normalized 2D gauss kernel array

    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)

    x, y = np.mgrid[-size:size + 1, -sizey:sizey + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(sizey)))

    return g / g.sum()


def interp_array(sci_data, dq_data, n_size):
    """
    Short Summary
    -------------
    For bad DQ values as given by the input DQ array, do a
    bilinear interpolation over the corresponding SCI values and return
    the interpolated SCI array.

    Parameters
    ----------
    sci_data: 2D float array
        original SCI image to interpolate over
    dq_data: 2D int array
        corresponding DQ image
    n_size: int
        size of the interpolation box

    Returns
    -------
    sci_data: 2D float array
        interpolated SCI image
    """

    wh_bad_dq = np.where(np.bitwise_and(dq_data, DO_NOT_USE))
    num_bad_dq = len(wh_bad_dq[0])

    # Create array of locations of bad DQ pixels
    bad_dq = np.zeros((num_bad_dq, 2), dtype=np.int16)
    for ii in np.arange(num_bad_dq):
        bad_dq[ii] = wh_bad_dq[0][ii], wh_bad_dq[1][ii]

    # Loop over bad pixels, filling in missing values with interpolated values
    for jj in np.arange(num_bad_dq):
        ga = create_griddata_array(sci_data, bad_dq[jj], n_size)

        x = bad_dq[jj][1]
        y = bad_dq[jj][0]

        # Linearly interpolate using scipy's griddata to fill in missing value
        sci_data[y, x] = griddata(ga[:, 0:2], ga[:, 2], [(y, x)], method='linear')

        # For those interpolations just done that result in a nan (because
        #    there may be too few pixels), check and redo with 'nearest'
        if np.isnan(sci_data[y, x]):
            sci_data[y, x] = griddata(ga[:, 0:2], ga[:, 2], [(y, x)], method='nearest')

    return sci_data


def create_griddata_array(sci_data, pixel, n_size):
    """
    Short Summary
    -------------
    Create interpolation array for input to scipy's griddata. This array
    consists of the coordinates and the pixel value for each of
    pixels neighboring the input pixel.

    Parameters
    ----------
    sci_data: 2D float array
        original SCI image

    pixel: int, int
        y, x coordinates of pixel to interpolate over

    n_size: int
        size of the interpolation box

    Returns
    -------
    interp_arr: int, int, float
        pixel coords, pixel value for each pixel neighboring the input pixel

    """
    xdim = sci_data.shape[1]
    ydim = sci_data.shape[0]

    # Generate neighborhood limits
    xmin = max(0, pixel[1] - n_size)
    ymin = max(0, pixel[0] - n_size)
    xmax = min(xdim - n_size, pixel[1] + n_size)
    ymax = min(ydim - n_size, pixel[0] + n_size)

    # Make a list for neighboring pixels, containing:
    # 1. coordinates for up to (2*N_SIZE+1)^2-1 neighbors, accounting for edges
    # 2. SCI data

    interp_list = []
    for x in range(xmin, xmax + 1):
        for y in range(ymin, ymax + 1):
            interp_list.append([y, x, sci_data[y, x]])

    # Remove identity element (central pixel)
    try:
        interp_list.remove([pixel[0], pixel[1], sci_data[pixel[0], pixel[1]]])
    except ValueError:
        pass

    interp_arr = np.asarray(interp_list)  # griddata requires an array

    return interp_arr


def get_index_range(offset, length):
    """
    Short Summary
    -------------
    Get the initial and final indices for the given offset and array length:
    For offset <= 0: i = 0,  f = length - abs(offset)
    For offset > 0: i = offset,  f = length

    Parameters
    ----------
    offset: integer
        offset

    length: integer
        length of (1D) array

    Returns
    -------
        initial index : integer
        final index : integer
    """
    if offset > 0:
        initial_1 = 0
        final_1 = length - offset
        initial_2 = offset
        final_2 = length
    else:
        initial_1 = abs(offset)
        final_1 = length
        initial_2 = 0
        final_2 = length - abs(offset)

    return initial_1, final_1, initial_2, final_2


def get_overlap(sci_int_1, sci_int_2, nom_off_x, nom_off_y):
    """
    Short Summary
    -------------
    Apply nominal offsets (of image #2 relative to image #1) to determine
    the overlap in interpolated images.

    Long Summary
    -------------
    Apply nominal offsets (of image #2 relative to image #1) to determine
    the overlap in interpolated images.  The resulting two subarrays are
    the pixels common to both. In other words, image #2 is shifted onto
    the frame of image #1, with the dimensions of the subarrays equal to
    the dimensions of the overlap.

    To illustrate with pseudocode for a 1D array with length 'length': for
    a given offset 'off', the resulting initial and final indices, and
    the elements of arrays indexed are:

    for off < 0 : ix = 0 and final_x = length - abs(off)
       subarray indices: sub_1[0: length - abs(off)]
       subarray indices: sub_2[0: length - abs(off)]

    for off = 0 : ix = 0 ; final_x = length)
       subarray indices: sub_1[0: length]
       subarray indices: sub_2[0: length]

    for off > 0 : ix = off ; final_x = length
       subarray indices: sub_1[0: length - off]
       subarray indices: sub_2[off: length]

    Parameters
    ----------
    sci_int_1: 2d float array
        interpolated SCI array for image 1
    sci_int_2: 2d float array
        interpolated SCI array for image 2
    nom_off_x: integer
        nominal offset in x-direction
    nom_off_y: integer
        nominal offset in y-direction

    Returns
    -------
    sub_1: 2d float array
        overlapping subarray for interpolated image 1
    sub_2: 2d float array
        overlapping subarray for interpolated image 2
    """

    # From the nominal offsets, determine array indices to shift image #2
    #     onto frame #1
    initial_x_1, final_x_1, initial_x_2, final_x_2 = get_index_range(nom_off_x, sci_int_2.shape[1])
    initial_y_1, final_y_1, initial_y_2, final_y_2 = get_index_range(nom_off_y, sci_int_2.shape[0])

    sub_1 = sci_int_1[initial_y_1:final_y_1, initial_x_1:final_x_1]
    sub_2 = sci_int_2[initial_y_2:final_y_2, initial_x_2:final_x_2]

    return sub_1, sub_2


def calc_refined_offsets(sci_nai_1, sci_nai_2, off_x, off_y, psf_size):
    """
    Short Summary
    -------------
    Get the overlap of the 2 images (based on the offsets), and
    calculate the two dimensional cross correlation image between 2 image subarrays.
    Then we slice on the a subarray around the peak of the cross correlation image and
    find the first moment. The first moment provides a high S/N measurement of the offset
    between the two images.


    Parameters
    ----------
    sci_nai_1: 2d float array
        nominally aligned, interpolated SCI subarray for image 1
    sci_nai_2: 2d float array
        nominally aligned, interpolated SCI subarray for image 2
    off_x: integer
        offset in x-direction
    off_y: integer
        offset in y-direction
    psf_size: integer
        The worst case PSF size

    Returns
    -------
    refined_x: float
        The refined value of the x offset
    refined_y: float
        The refined value of the y offset

    """
    centroid_size = 3
    sub_1, sub_2 = get_overlap(sci_nai_1, sci_nai_2, off_x, off_y)
    num_pix = sub_1.shape[0] * sub_1.shape[1]  # Number of overlapping pixels

    # Raise (fatal) exception if there are no overlapping pixels
    if num_pix == 0:
        log.error('Applying offsets to image #2 results in 0 overlapping pix')
        raise RuntimeWarning('No overlapping pixels in 2 images')

    # Set limits for subarrays, centered on the overlap and +/- psf half width,
    #   taking edges into account
    xlen = sub_1.shape[1]
    ylen = sub_1.shape[0]
    xcen = xlen / 2 + 1
    ycen = ylen / 2 + 1

    xmin = int(max(0, xcen - psf_size))
    xmax = int(min(xlen - 1, xcen + psf_size))
    ymin = int(max(0, ycen - psf_size))
    ymax = int(min(ylen - 1, ycen + psf_size))

    # Create subarrays using these limits
    sub_1_sub = sub_1[ymin:ymax, xmin:xmax]
    sub_2_sub = sub_2[ymin:ymax, xmin:xmax]
    # Create the cross correlation image
    cross_cor = correlate2d(sub_2_sub - gaussian_filter(sub_2_sub, 5), sub_1_sub -
                            gaussian_filter(sub_1_sub, 5))
    maximum_pixel = np.unravel_index(np.argmax(cross_cor), cross_cor.shape)

    ymax = maximum_pixel[0] - sub_1_sub.shape[0] + 1
    xmax = maximum_pixel[1] - sub_1_sub.shape[1] + 1
    # Slice out a box center on the peak of the cross correlation image. The centroid of this box will give a
    # accurate estimate of the x and y offsets.
    central_cutout = cross_cor[maximum_pixel[0] - centroid_size:maximum_pixel[0] + centroid_size + 1,
                               maximum_pixel[1] - centroid_size:maximum_pixel[1] + centroid_size + 1]
    centroid = center_of_mass(central_cutout)
    refined_x = xmax + centroid[1] - centroid_size
    refined_y = ymax + centroid[0] - centroid_size
    return refined_x, refined_y
