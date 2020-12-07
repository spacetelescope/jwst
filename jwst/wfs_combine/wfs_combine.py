import logging

import numpy as np
from scipy.interpolate import griddata
from scipy.signal import convolve
from scipy import mgrid

from .. import datamodels
from ..datamodels import dqflags


log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# The following should probably be in a ref file instead
OFF_DELTA = 2  # for searching +/-2 around nominal offsets
PSF_SIZE = 200  # half width of psf for 12 waves
BLUR_SIZE = 10  # size of gaussian kernel for convolution
N_SIZE = 2  # size of neighborhood in create_griddata_array


DO_NOT_USE = dqflags.pixel["DO_NOT_USE"]


class DataSet:
    """
    Two dithered input wavefront sensing images to be combined

    """

    def __init__(self, infile_1, infile_2, outfile, do_refine):
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

        if (self.input_1.data.shape != self.input_2.data.shape):
            log.error('Incompatible sizes for input files')

        log.info('File 1 to combine: %s', infile_1)
        log.info('File 2 to combine: %s', infile_2)
        log.info('Output file: %s', outfile)
        log.info('do_refine: %s', do_refine)

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

        # Input SCI arrays may have nan's so replace with 0's to prevent
        # later annoyances (hopefully this can be removed later)
        self.input_1.data[np.isnan(self.input_1.data)] = 0.
        self.input_2.data[np.isnan(self.input_2.data)] = 0.

        im_1_a = self.input_1.copy()  # Aligned image #1 (already aligned)
        im_2_a = self.create_aligned_2()  # Aligned image #2

        # Create and populate extensions for combined data
        data_c, dq_c, err_c = self.create_combined(im_1_a, im_2_a)

        # Create a new model using the combined arrays...
        new_model = datamodels.ImageModel(data=data_c, dq=dq_c, err=err_c)
        new_model.update(self.input_1)

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
           centroid +/- PSF_SIZE and adding the BLUR_SIZE, taking the edges
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

        if self.do_refine:
            # 1. Create smoothed image of input SCI data of image #1
            # 1a. create image to smooth by first setting bad DQ pixels equal
            #     to mean of good pixels
            data_1 = self.input_1.data.astype(np.float)
            bad1 = np.bitwise_and(self.input_1.dq, DO_NOT_USE).astype(np.bool)
            data_1[bad1] = data_1[~bad1].mean()

            # 1b. Create smoothed image by smoothing this 'repaired' image
            g = gauss_kern(BLUR_SIZE, sizey=None)
            s_data_1 = convolve(data_1, g, mode='valid')

            # 2. Find approximate center of PSF in umsmoothed frame by taking
            #    all pixels in smoothed image exceeding 50% of the maximum
            #    of the smoothed image, and taking the mean of the coordinates
            #    of these pixels. Add BLUR_SIZE to take smoothing into account
            wh_data_hi = np.where(s_data_1 > 0.5 * s_data_1.max())

            ctrd_x = wh_data_hi[1].mean() + BLUR_SIZE
            ctrd_y = wh_data_hi[0].mean() + BLUR_SIZE

            log.info('Approximate centroid of image 1 PSF has x,y : %s %s',
                     round(ctrd_x), round(ctrd_y))

            # 3. Set limits of the subarrays (in frames of input data)
            #    for interpolation by taking this centroid +/- PSF_SIZE
            #    and adding BLUR_SIZE, taking edges into account
            xmin = int(round(max(0, ctrd_x - PSF_SIZE)))
            ymin = int(round(max(0, ctrd_y - PSF_SIZE)))
            xmax = int(round(min(self.input_1.data.shape[1], ctrd_x + PSF_SIZE)))
            ymax = int(round(min(self.input_1.data.shape[0], ctrd_y + PSF_SIZE)))

            # 3a. Set subarrays and interpolate over bad pixels
            data_sub_1 = self.input_1.data[ymin: ymax, xmin: xmax]
            dq_sub_1 = self.input_1.dq[ymin: ymax, xmin: xmax]
            sci_int_1 = interp_array(data_sub_1, dq_sub_1)

            data_sub_2 = self.input_2.data[ymin: ymax, xmin: xmax]
            dq_sub_2 = self.input_2.dq[ymin: ymax, xmin: xmax]
            sci_int_2 = interp_array(data_sub_2, dq_sub_2)

            # 4. Determine overlap of these interpolated images, and
            #    return nominally aligned, interpolated images
            sci_nai_1, sci_nai_2 = get_overlap(sci_int_1, sci_int_2,
                                               self.off_x, self.off_y)

            # 5. Around this nominal alignment, get refined (delta) offsets
            ref_del_off_x, ref_del_off_y = optimize_offs(sci_nai_1, sci_nai_2)

            log.info('From the refined offsets calculation,'
                     'the x,y changes in ofsets are: %s %s',
                     ref_del_off_x, ref_del_off_y)

            # 6. Add the refined delta offsets to the nominal offsets
            self.off_x += ref_del_off_x
            self.off_y += ref_del_off_y

            log.info('Values for the refined offsets are, for x,y : %s %s',
                     self.off_x, self.off_y)

        log.info(f"Final x,y offset in pixels: {self.off_x} {self.off_y}")

        # Do final alignment for original (not interpolated) image #2
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
        xcen = int(self.input_1.data.shape[1]/2)
        ycen = int(self.input_1.data.shape[0]/2)

        radec = tr1(xcen, ycen)
        pixels = tr2(radec[0], radec[1])
        off_x = pixels[0] - xcen
        off_y = pixels[1] - ycen

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

        data1 = image1.data.astype(np.float)
        data2 = image2.data.astype(np.float)
        dq1 = image1.dq.copy()
        dq2 = image2.dq.copy()
        err1 = image1.err.copy()
        err2 = image2.err.copy()

        # Create boolean arrays of bad pixels in each input image
        bad1 = np.bitwise_and(dq1, DO_NOT_USE).astype(np.bool)
        good1 = ~bad1
        bad2 = np.bitwise_and(dq2, DO_NOT_USE).astype(np.bool)
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

        return data_comb, dq_comb, err_comb

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

        ai_x, af_x = get_index_range(self.off_x, a.shape[1])
        ai_y, af_y = get_index_range(self.off_y, a.shape[0])

        bi_x = a.shape[1] - af_x  # For output, x-direction's initial channel
        bf_x = a.shape[1] - ai_x  # ...and final channel

        bi_y = a.shape[0] - af_y  # For output, y-direction's initial channel
        bf_y = a.shape[0] - ai_y  # ...and final channel

        b = np.zeros(a.shape)
        b[bi_y:bf_y, bi_x:bf_x] = a[ai_y:af_y, ai_x:af_x]

        return b


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

    x, y = mgrid[-size:size + 1, -sizey:sizey + 1]
    g = np.exp(-(x**2 / float(size) + y**2 / float(sizey)))

    return g / g.sum()


def interp_array(sci_data, dq_data):
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
        ga = create_griddata_array(sci_data, bad_dq[jj])

        x = bad_dq[jj][1]
        y = bad_dq[jj][0]

        # Linearly interpolate using scipy's griddata to fill in missing value
        sci_data[y, x] = griddata(ga[:, 0:2], ga[:, 2], [(y, x)], method='linear')

        # For those interpolations just done that result in a nan (because
        #    there may be too few pixels), check and redo with 'nearest'
        if np.isnan(sci_data[y, x]):
            sci_data[y, x] = griddata(ga[:, 0:2], ga[:, 2], [(y, x)], method='nearest')

    return sci_data


def create_griddata_array(sci_data, pixel):
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

    Returns
    -------
    interp_arr: int, int, float
        pixel coords, pixel value for each pixel neighboring the input pixel

    """
    xdim = sci_data.shape[1]
    ydim = sci_data.shape[0]

    # Generate neighborhood limits
    xmin = max(0, pixel[1] - N_SIZE)
    ymin = max(0, pixel[0] - N_SIZE)
    xmax = min(xdim - N_SIZE, pixel[1] + N_SIZE)
    ymax = min(ydim - N_SIZE, pixel[0] + N_SIZE)

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
    i: integer
        initial index
    f: integer
        final index
    """

    i = int((abs(offset) + offset) / 2)
    f = length + int((-abs(offset) + offset) / 2)

    return i, f


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

    for off < 0 : ix = 0 and fx = length - abs(off)
       subarray indices: sub_1[0: length - abs(off)]
       subarray indices: sub_2[0: length - abs(off)]

    for off = 0 : ix = 0 ; fx = length)
       subarray indices: sub_1[0: length]
       subarray indices: sub_2[0: length]

    for off > 0 : ix = off ; fx = length
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
    ix, fx = get_index_range(nom_off_x, sci_int_2.shape[1])
    iy, fy = get_index_range(nom_off_y, sci_int_2.shape[0])

    sub_1 = sci_int_1[:fy - iy, :fx - ix]
    sub_2 = sci_int_2[iy:fy, ix:fx]

    return sub_1, sub_2


def optimize_offs(sci_nai_1, sci_nai_2):
    """
    Short Summary
    -------------
    Calculate refined values for the offsets. Find the nearby integer
    relative offsets between the nominally aligned interpolated images
    which minimizes the square of the difference between the 2 images.

    Parameters
    ----------
    sci_nai_1: 2d float array
        nominally aligned, interpolated SCI array for image 1
    sci_nai_2: 2d float array
        nominally aligned, interpolated SCI array for image 2

    Returns
    -------
    max_off_x: integer
        offset in x-direction which minimizes array difference
    max_off_y: integer
        offset in y-direction which minimizes array difference

    """

    max_diff = 0.

    for off_y in range(-OFF_DELTA, OFF_DELTA + 1):
        for off_x in range(-OFF_DELTA, OFF_DELTA + 1):
            this_diff = calc_cor_coef(sci_nai_1, sci_nai_2, off_x, off_y)

            if (this_diff > max_diff):
                max_diff = this_diff
                max_off_x = off_x
                max_off_y = off_y

    return max_off_x, max_off_y


def calc_cor_coef(sci_nai_1, sci_nai_2, off_x, off_y):
    """
    Short Summary
    -------------
    Get the overlap of the 2 images (based on the offsets), and
    calculate the correlation coefficient between 2 image subarrays.
    (We may want something more sophisticated later)

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

    Returns
    -------
    mean_diff: float
        mean of correlation coefficient array

    """

    sub_1, sub_2 = get_overlap(sci_nai_1, sci_nai_2, off_x, off_y)
    num_pix = sub_1.shape[0] * sub_1.shape[1]  # Number of overlapping pixels

    # Raise (fatal) exception if there are no overlapping pixels
    if (num_pix == 0):
        log.error('Applying offsets to image #2 results in 0 overlapping pix')
        raise RuntimeWarning('No overlapping pixels in 2 images')

    # Set limits for subarrays, centered on the overlap and +/- psf half width,
    #   taking edges into account
    xlen = sub_1.shape[1]
    ylen = sub_1.shape[0]
    xcen = xlen / 2 + 1
    ycen = ylen / 2 + 1

    xmin = int(max(0, xcen - PSF_SIZE))
    xmax = int(min(xlen - 1, xcen + PSF_SIZE))
    ymin = int(max(0, ycen - PSF_SIZE))
    ymax = int(min(ylen - 1, ycen + PSF_SIZE))

    # Create subarrays using these limits, and make them
    #   1D for numpy's correlation coefficient
    sub_1_sub = sub_1[ymin:ymax, xmin:xmax]
    sub_2_sub = sub_2[ymin:ymax, xmin:xmax]
    sub_1_sub = np.ravel(sub_1_sub)
    sub_2_sub = np.ravel(sub_2_sub)

    cor_coef = np.corrcoef(sub_1_sub, sub_2_sub)

    return cor_coef.mean()
