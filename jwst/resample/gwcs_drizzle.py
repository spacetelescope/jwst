import numpy as np

from drizzle import util
from drizzle import doblot
from drizzle import cdrizzle
from . import resample_utils

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

class GWCSDrizzle:
    """
    Combine images using the drizzle algorithm
    """
    def __init__(self, product, outwcs=None, single=False,
                 wt_scl="exptime", pixfrac=1.0, kernel="square",
                 fillval="INDEF"):
        """
        Create a new Drizzle output object and set the drizzle parameters.

        Parameters
        ----------

        product : str, optional
            A data model containing results from a previous run. The three
            extensions SCI, WHT, and CTX contain the combined image, total counts
            and image id bitmap, repectively. The WCS of the combined image is
            also read from the SCI extension.

        outwcs : `gwcs.WCS`
            The world coordinate system (WCS) of the resampled image.  If not
            provided, the WCS is taken from product.

        wt_scl : str, optional
            How each input image should be scaled. The choices are `exptime`
            which scales each image by its exposure time, `expsq` which scales
            each image by the exposure time squared, or an empty string, which
            allows each input image to be scaled individually.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel : str, optional
            The name of the kernel used to combine the inputs. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: "square", "gaussian", "point", "tophat", "turbo", "lanczos2",
            and "lanczos3". The square kernel is the default.

        fillval : str, otional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of INDEF does not set a value.
        """

        # Initialize the object fields
        self.outsci = None
        self.outwht = None
        self.outcon = None

        self.outexptime = 0.0
        self.uniqid = 0

        self.wt_scl = wt_scl
        self.kernel = kernel
        self.fillval = fillval
        self.pixfrac = pixfrac

        self.sciext = "SCI"
        self.whtext = "WHT"
        self.conext = "CON"

        out_units = "cps"

        self.outexptime = product.meta.resample.product_exposure_time or 0.0

        self.outsci = product.data
        if outwcs:
            self.outwcs = outwcs
        else:
            self.outwcs = product.meta.wcs

        self.outwht = product.wht
        self.outcon = product.con

        if self.outcon.ndim == 2:
            self.outcon = np.reshape(self.outcon, (1,
                                     self.outcon.shape[0],
                                     self.outcon.shape[1]))

        elif self.outcon.ndim == 3:
            pass

        else:
            raise ValueError("Drizzle context image has wrong dimensions: \
                {0}".format(product))

        # Check field values
        if not self.outwcs:
            raise ValueError("Either an existing file or wcs must be supplied")

        if util.is_blank(self.wt_scl):
            self.wt_scl = ''
        elif self.wt_scl != "exptime" and self.wt_scl != "expsq":
            raise ValueError("Illegal value for wt_scl: %s" % self.wt_scl)

        if out_units == "counts":
            np.divide(self.outsci, self.outexptime, self.outsci)
        elif out_units != "cps":
            raise ValueError("Illegal value for out_units: %s" % out_units)

    def add_image(self, insci, inwcs, inwht=None, xmin=0, xmax=0, ymin=0, ymax=0,
                  expin=1.0, in_units="cps", wt_scl=1.0):
        """
        Combine an input image with the output drizzled image.

        Instead of reading the parameters from a fits file, you can set
        them by calling this lower level method. `Add_fits_file` calls
        this method after doing its setup.

        Parameters
        ----------

        insci : array
            A 2d numpy array containing the input image to be drizzled.
            it is an error to not supply an image.

        inwcs : wcs
            The world coordinate system of the input image. This is
            used to convert the pixels to the output coordinate system.

        inwht : array, optional
            A 2d numpy array containing the pixel by pixel weighting.
            Must have the same dimenstions as insci. If none is supplied,
            the weghting is set to one.

        xmin : float, optional
            This and the following three parameters set a bounding rectangle
            on the output image. Only pixels on the output image inside this
            rectangle will have their flux updated. Xmin sets the minimum value
            of the x dimension. The x dimension is the dimension that varies
            quickest on the image. If the value is zero or less, no minimum will
            be set in the x dimension. All four parameters are zero based,
            counting starts at zero.

        xmax : float, optional
            Sets the maximum value of the x dimension on the bounding box
            of the ouput image. If the value is zero or less, no maximum will
            be set in the x dimension.

        ymin : float, optional
            Sets the minimum value in the y dimension on the bounding box. The
            y dimension varies less rapidly than the x and represents the line
            index on the output image. If the value is zero or less, no minimum
            will be set in the y dimension.

        ymax : float, optional
            Sets the maximum value in the y dimension. If the value is zero or
            less, no maximum will be set in the y dimension.

        expin : float, optional
            The exposure time of the input image, a positive number. The
            exposure time is used to scale the image if the units are counts and
            to scale the image weighting if the drizzle was initialized with
            wt_scl equal to "exptime" or "expsq."

        in_units : str, optional
            The units of the input image. The units can either be "counts"
            or "cps" (counts per second.) If the value is counts, before using
            the input image it is scaled by dividing it by the exposure time.

        wt_scl : float, optional
            If drizzle was initialized with wt_scl left blank, this value will
            set a scaling factor for the pixel weighting. If drizzle was
            initialized with wt_scl set to "exptime" or "expsq", the exposure time
            will be used to set the weight scaling and the value of this parameter
            will be ignored.
        """
        insci = insci.astype(np.float32)

        if inwht is None:
            inwht = np.ones(insci.shape, dtype=insci.dtype)
        else:
            inwht = inwht.astype(np.float32)

        if self.wt_scl == "exptime":
            wt_scl = expin
        elif self.wt_scl == "expsq":
            wt_scl = expin * expin

        wt_scl = 1.0  # hard-coded for JWST count-rate data
        self.increment_id()

        dodrizzle(insci, inwcs, inwht, self.outwcs, self.outsci, self.outwht,
                  self.outcon, expin, in_units, wt_scl, uniqid=self.uniqid,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                  pixfrac=self.pixfrac, kernel=self.kernel, fillval=self.fillval)

    def blot_image(self, blotwcs, interp='poly5', sinscl=1.0):
        """
        Resample the output image using an input world coordinate system.

        Parameters
        ----------

        blotwcs : wcs
            The world coordinate system to resample on.

        interp : str, optional
            The type of interpolation used in the resampling. The
            possible values are "nearest" (nearest neighbor interpolation),
            "linear" (bilinear interpolation), "poly3" (cubic polynomial
            interpolation), "poly5" (quintic polynomial interpolation),
            "sinc" (sinc interpolation), "lan3" (3rd order Lanczos
            interpolation), and "lan5" (5th order Lanczos interpolation).

        sincscl : float, optional
            The scaling factor for sinc interpolation.
        """

        util.set_pscale(blotwcs)
        self.outsci = doblot.doblot(self.outsci, self.outwcs, blotwcs,
                                    1.0, interp=interp, sinscl=sinscl)

        self.outwcs = blotwcs

    def increment_id(self):
        """
        Increment the id count and add a plane to the context image if needed

        Drizzle tracks which input images contribute to the output image
        by setting a bit in the corresponding pixel in the context image.
        The uniqid indicates which bit. So it must be incremented each time
        a new image is added. Each plane in the context image can hold 32 bits,
        so after each 32 images, a new plane is added to the context.
        """

        # Compute what plane of the context image this input would
        # correspond to:
        planeid = int(self.uniqid / 32)

        # Add a new plane to the context image if planeid overflows

        if self.outcon.shape[0] == planeid:
            plane = np.zeros_like(self.outcon[0])
            plane = plane.reshape((1, plane.shape[0], plane.shape[1]))
            self.outcon = np.concatenate((self.outcon, plane))

        # Increment the id
        self.uniqid += 1


def dodrizzle(insci, input_wcs, inwht, output_wcs, outsci, outwht, outcon,
              expin, in_units, wt_scl, uniqid=1, xmin=0, xmax=0, ymin=0, ymax=0,
              pixfrac=1.0, kernel='square', fillval="INDEF"):
    """
    Low level routine for performing 'drizzle' operation on one image.

    The interface is compatible with STScI code. All images are Python
    ndarrays, instead of filenames. File handling (input and output) is
    performed by the calling routine.

    Parameters
    ----------

    insci : 2d array
        A 2d numpy array containing the input image to be drizzled.

    input_wcs : gwcs.WCS object
        The world coordinate system of the input image.

    inwht : 2d array
        A 2d numpy array containing the pixel by pixel weighting.
        Must have the same dimensions as insci. If none is supplied,
        the weghting is set to one.

    output_wcs : gwcs.WCS object
        The world coordinate system of the output image.

    outsci : 2d array
        A 2d numpy array containing the output image produced by
        drizzling. On the first call it should be set to zero.
        Subsequent calls it will hold the intermediate results

    outwht : 2d array
        A 2d numpy array containing the output counts. On the first
        call it should be set to zero. On subsequent calls it will
        hold the intermediate results.

    outcon : 2d or 3d array, optional
        A 2d or 3d numpy array holding a bitmap of which image was an input
        for each output pixel. Should be integer zero on first call.
        Subsequent calls hold intermediate results.

    expin : float
        The exposure time of the input image, a positive number. The
        exposure time is used to scale the image if the units are counts.

    in_units : str
        The units of the input image. The units can either be "counts"
        or "cps" (counts per second.)

    wt_scl : float
        A scaling factor applied to the pixel by pixel weighting.

    uniqid : int, optional
        The id number of the input image. Should be one the first time
        this function is called and incremented by one on each subsequent
        call.

    xmin : float, optional
        This and the following three parameters set a bounding rectangle
        on the input image. Only pixels on the input image inside this
        rectangle will have their flux added to the output image. Xmin
        sets the minimum value of the x dimension. The x dimension is the
        dimension that varies quickest on the image. If the value is zero,
        no minimum will be set in the x dimension. All four parameters are
        zero based, counting starts at zero.

    xmax : float, optional
        Sets the maximum value of the x dimension on the bounding box
        of the input image. If the value is zero, no maximum will
        be set in the x dimension, the full x dimension of the output
        image is the bounding box.

    ymin : float, optional
        Sets the minimum value in the y dimension on the bounding box. The
        y dimension varies less rapidly than the x and represents the line
        index on the input image. If the value is zero, no minimum  will be
        set in the y dimension.

    ymax : float, optional
        Sets the maximum value in the y dimension. If the value is zero, no
        maximum will be set in the y dimension,  the full x dimension
        of the output image is the bounding box.

    pixfrac : float, optional
        The fraction of a pixel that the pixel flux is confined to. The
        default value of 1 has the pixel flux evenly spread across the image.
        A value of 0.5 confines it to half a pixel in the linear dimension,
        so the flux is confined to a quarter of the pixel area when the square
        kernel is used.

    kernel: str, optional
        The name of the kernel used to combine the input. The choice of
        kernel controls the distribution of flux over the kernel. The kernel
        names are: "square", "gaussian", "point", "tophat", "turbo", "lanczos2",
        and "lanczos3". The square kernel is the default.

    fillval: str, optional
        The value a pixel is set to in the output if the input image does
        not overlap it. The default value of INDEF does not set a value.

    Returns
    -------
    A tuple with three values: a version string, the number of pixels
    on the input image that do not overlap the output image, and the
    number of complete lines on the input image that do not overlap the
    output input image.

    """

    # Insure that the fillval parameter gets properly interpreted for use with tdriz
    if util.is_blank(str(fillval)):
        fillval = 'INDEF'
    else:
        fillval = str(fillval)

    if in_units == 'cps':
        expscale = 1.0
    else:
        expscale = expin

    if (insci.dtype > np.float32):
        insci = insci.astype(np.float32)

    # Add input weight image if it was not passed in
    if inwht is None:
        inwht = np.ones_like(insci)

    if xmax is None or xmax == xmin:
        xmax = insci.shape[1]
    if ymax is None or ymax == ymin:
        ymax = insci.shape[0]

    # Compute what plane of the context image this input would
    # correspond to:
    planeid = int((uniqid - 1) / 32)

    # Check if the context image has this many planes
    if outcon.ndim == 3:
        nplanes = outcon.shape[0]
    elif outcon.ndim == 2:
        nplanes = 1
    else:
        nplanes = 0

    if nplanes <= planeid:
        raise IndexError("Not enough planes in drizzle context image")

    # Alias context image to the requested plane if 3d
    if outcon.ndim == 3:
        outcon = outcon[planeid]

    # Compute the mapping between the input and output pixel coordinates
    # for use in drizzle.cdrizzle.tdriz
    pixmap = resample_utils.calc_gwcs_pixmap(input_wcs, output_wcs, insci.shape)
    # inwht[np.isnan(pixmap[:,:,0])] = 0.

    log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
    log.debug(f"Input Sci shape: {insci.shape}")
    log.debug(f"Output Sci shape: {outsci.shape}")

    # y_mid = pixmap.shape[0] // 2
    # x_mid = pixmap.shape[1] // 2
    # print("x slice: ", pixmap[y_mid,:,0])
    # print("y slice: ", pixmap[:,x_mid,1])
    # print("insci: ", insci)
    # Call 'drizzle' to perform image combination

    log.info(f"Drizzling {insci.shape} --> {outsci.shape}")

    _vers, nmiss, nskip = cdrizzle.tdriz(
        insci, inwht, pixmap,
        outsci, outwht, outcon,
        uniqid=uniqid,
        xmin=xmin, xmax=xmax,
        ymin=ymin, ymax=ymax,
        pixfrac=pixfrac,
        kernel=kernel,
        in_units=in_units,
        expscale=expscale,
        wtscale=wt_scl,
        fillstr=fillval
        )
    return _vers, nmiss, nskip
