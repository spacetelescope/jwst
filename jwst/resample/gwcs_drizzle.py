
from drizzle.resample import Drizzle
from . import resample_utils

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class GWCSDrizzle(Drizzle):
    """
    Combine images using the drizzle algorithm

    """
    def __init__(self, output_wcs=None, pixfrac=1.0, kernel="square",
                 fillval="NAN", disable_ctx=False, n_max_images=None):
        """
        Create a new Drizzle output object and set the drizzle parameters.

        Parameters
        ----------

        product : DataModel
            A data model containing results from a previous run. The three
            extensions SCI, WHT, and CTX contain the combined image, total counts
            and image id bitmap, respectively. The WCS of the combined image is
            also read from the SCI extension.

        outwcs : `gwcs.WCS`
            The world coordinate system (WCS) of the resampled image.  If not
            provided, the WCS is taken from product.

        pixfrac : float, optional
            The fraction of a pixel that the pixel flux is confined to. The
            default value of 1 has the pixel flux evenly spread across the image.
            A value of 0.5 confines it to half a pixel in the linear dimension,
            so the flux is confined to a quarter of the pixel area when the square
            kernel is used.

        kernel : str, optional
            The name of the kernel used to combine the inputs. The choice of
            kernel controls the distribution of flux over the kernel. The kernel
            names are: "square", "gaussian", "point", "turbo", "lanczos2",
            and "lanczos3". The square kernel is the default.

        fillval : str, optional
            The value a pixel is set to in the output if the input image does
            not overlap it. The default value of NAN sets NaN values.
        """
        self.output_wcs = output_wcs
        self.pixfrac = pixfrac

        super().__init__(
            kernel=kernel,
            fillval=fillval,
            out_shape=output_wcs.array_shape,
            out_img=None,
            out_wht=None,
            out_ctx=None,
            exptime=0.0,
            begin_ctx_id=0,
            max_ctx_id=n_max_images,
            disable_ctx=disable_ctx
        )

    def add_image(self, insci, inwcs, pixmap=None, inwht=None,
                  xmin=0, xmax=0, ymin=0, ymax=0,
                  expin=1.0, in_units="cps", iscale=1.0):
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
            Must have the same dimensions as insci. If none is supplied,
            the weighting is set to one.

        pixmap : array, optional
            A 3D array that maps input image coordinates onto the output
            array. If provided, it can improve performance.

        xmin : int, optional
            This and the following three parameters set a bounding rectangle
            on the input image. Only pixels on the input image inside this
            rectangle will have their flux added to the output image. Xmin
            sets the minimum value of the x dimension. The x dimension is the
            dimension that varies quickest on the image. All four parameters
            are zero based, counting starts at zero.

        xmax : int, optional
            Sets the maximum value of the x dimension on the bounding box
            of the input image. If ``xmax = 0``, no maximum will
            be set in the x dimension (all pixels in a row of the input image
            will be resampled).

        ymin : int, optional
            Sets the minimum value in the y dimension on the bounding box. The
            y dimension varies less rapidly than the x and represents the line
            index on the input image.

        ymax : int, optional
            Sets the maximum value in the y dimension. If ``ymax = 0``,
            no maximum will be set in the y dimension (all pixels in a column
            of the input image will be resampled).

        expin : float, optional
            The exposure time of the input image, a positive number. The
            exposure time is used to scale the image if the units are counts and
            to scale the image weighting if the drizzle was initialized with
            wt_scl equal to "exptime" or "expsq."

        in_units : str, optional
            The units of the input image. The units can either be "counts"
            or "cps" (counts per second.) If the value is counts, before using
            the input image it is scaled by dividing it by the exposure time.

        iscale : float, optional
            A scale factor to be applied to pixel intensities of the
            input image before resampling.

        """
        # Compute the mapping between the input and output pixel coordinates
        # for use in drizzle.cdrizzle.tdriz
        if pixmap is None:
            pixmap = resample_utils.calc_gwcs_pixmap(
                inwcs,
                self.output_wcs,
                insci.shape
            )

        log.debug(f"Pixmap shape: {pixmap[:,:,0].shape}")
        log.debug(f"Input Sci shape: {insci.shape}")
        log.debug(f"Output Sci shape: {self.out_img.shape}")

        # Call 'drizzle' to perform image combination
        log.info(f"Drizzling {insci.shape} --> {self.out_img.shape}")

        super().add_image(
            data=insci,
            exptime=expin,
            pixmap=pixmap,
            scale=iscale,
            weight_map=inwht,
            wht_scale=1.0,  # hard-coded for JWST count-rate data
            pixfrac=self.pixfrac,
            in_units=in_units,
            xmin=xmin,
            xmax=xmax,
            ymin=ymin,
            ymax=ymax,
        )
