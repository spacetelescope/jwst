import numpy as np

from scipy import sparse

from stdatamodels.jwst import datamodels

from jwst.lib.winclip import get_clipped_pixels
from .sens1d import create_1d_sens

import logging
import warnings

from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def background_subtract(
    data,
    box_size=None,
    filter_size=(3, 3),
    sigma=3.0,
    exclude_percentile=30.0,
):
    """
    Apply a simple astropy background subtraction.

    Parameters
    ----------
    data : np.ndarray
        2D array of pixel values
    box_size : tuple
        Size of box in pixels to use for background estimation.
        If not set, defaults to 1/5 of the image size.
    filter_size : tuple
        Size of filter to use for background estimation
    sigma : float
        Sigma threshold for background clipping
    exclude_percentile : float
        Percentage of masked pixels above which box is excluded from background estimation

    Returns
    -------
    data : np.ndarray
        2D array of pixel values with background subtracted

    Notes
    -----
    Improper background subtraction in input _i2d image leads to extra flux
    in the simulated dispersed image, and was one cause of flux scaling issues
    in a previous version.
    """
    if box_size is None:
        box_size = (int(data.shape[0] / 5), int(data.shape[1] / 5))
    sigma_clip = SigmaClip(sigma=sigma)
    bkg_estimator = MedianBackground()
    bkg = Background2D(
        data,
        box_size,
        filter_size=filter_size,
        sigma_clip=sigma_clip,
        bkg_estimator=bkg_estimator,
        exclude_percentile=exclude_percentile,
    )

    return data - bkg.background


def _select_ids(source_id, all_ids):
    """
    Select the source IDs to be processed based on the input ID parameter.

    Parameters
    ----------
    source_id : int or list-like
        ID(s) of source to process. If None, all sources processed.
    all_ids : np.ndarray
        Array of all source IDs in the segmentation map

    Returns
    -------
    selected_IDs : list
        List of selected source IDs
    """
    if source_id is None:
        log.info(f"Loading all {len(all_ids)} sources from segmentation map")
        return all_ids

    elif isinstance(source_id, int):
        log.info(f"Loading single source {source_id} from segmentation map")
        return [source_id]

    elif isinstance(source_id, list) or isinstance(source_id, np.ndarray):
        log.info(
            f"Loading {len(source_id)} of {len(all_ids)} selected sources from segmentation map"
        )
        return list(source_id)
    else:
        raise ValueError("ID must be an integer or a list of integers")


class Observation:
    """Define an observation leading to a single grism image."""

    def __init__(
        self,
        direct_image,
        segmap_model,
        grism_wcs,
        filter_name,
        source_id=None,
        extrapolate_sed=False,
        boundaries=None,
        offsets=None,
        max_cpu=1,
        oversample_factor=2,
    ):
        """
        Initialize all data and metadata for a given observation.

        Creates lists of
        direct image pixel values for selected objects.

        Parameters
        ----------
        direct_image : str
            File name containing direct imaging data
        segmap_model : `jwst.datamodels.ImageModel`
            Segmentation map model
        grism_wcs : gwcs object
            WCS object from grism image
        filter_name : str
            Filter name
        source_id : int, optional, default 0
            ID of source to process. If 0, all sources processed.
        extrapolate_sed : bool, optional, default False
            Flag indicating whether to extrapolate wavelength range of SED
        boundaries : list, optional, default []
            Start/Stop coordinates of the FOV within the larger seed image.
        offsets : list, optional, default [0,0]
            Offset values for x and y axes
        max_cpu : int, optional, default 1
            Max number of cpu's to use when multiprocessing
        """
        if boundaries is None:
            boundaries = []
        if offsets is None:
            offsets = [0, 0]
        # Load all the info for this grism mode
        self.seg_wcs = segmap_model.meta.wcs
        self.grism_wcs = grism_wcs
        self.dir_image_name = direct_image
        self.seg = segmap_model.data
        all_ids = np.array(list(set(np.ravel(self.seg))))
        self.source_ids = _select_ids(source_id, all_ids)
        self.filter = filter_name
        self.pivlam = float(self.filter[1:4]) / 100.0
        self.max_cpu = max_cpu
        self.oversample_factor = oversample_factor
        self.xoffset = offsets[0]
        self.yoffset = offsets[1]

        # make the direct image
        with datamodels.open(self.dir_image_name) as model:
            dimage = model.data
            self.dimage = background_subtract(dimage)

        # Set the limits of the dispersed image to be simulated
        if len(boundaries) == 0:
            log.debug("No boundaries passed.")
            self.xstart = 0
            self.xend = self.xstart + self.seg.shape[0] - 1
            self.ystart = 0
            self.yend = self.ystart + self.seg.shape[1] - 1
        else:
            self.xstart, self.xend, self.ystart, self.yend = boundaries
        self.dims = (self.yend - self.ystart + 1, self.xend - self.xstart + 1)
        log.debug(f"Using simulated image size of ({self.dims[1]}, {self.dims[0]}).")

        # Allow for SED extrapolation
        self.extrapolate_sed = extrapolate_sed
        if self.extrapolate_sed:
            log.warning("SED Extrapolation turned on.")

        # Create pixel lists for sources labeled in segmentation map
        self._create_pixel_list()

        # Initialize the list of slits
        self.simul_slits = datamodels.MultiSlitModel()
        self.simul_slits_order = []
        self.simul_slits_sid = []

    def _create_pixel_list(self):
        """Create flat lists of pixels to be dispersed, grouped per object ID."""
        self.xs = []
        self.ys = []
        self.source_ids_per_pixel = []
        self.fluxes = []
        for source_id in self.source_ids:
            ys, xs = np.nonzero(self.seg == source_id)
            self.xs.extend(xs)
            self.ys.extend(ys)
            self.source_ids_per_pixel.extend([source_id] * len(xs))
            self.fluxes.extend(self.dimage[ys, xs])
        self.xs = np.array(self.xs)
        self.ys = np.array(self.ys)
        self.fluxes = np.array(self.fluxes)
        self.source_ids_per_pixel = np.array(self.source_ids_per_pixel)

    def disperse_one_order(
        self,
        order,
        wmin,
        wmax,
        sens_waves,
        sens_resp,
    ):
        """
        Compute dispersed pixel values for all sources identified in the segmentation map.

        Parameters
        ----------
        order : int
            Spectral order number to process
        wmin : float
            Minimum wavelength for dispersed spectra
        wmax : float
            Maximum wavelength for dispersed spectra
        sens_waves : float array
            Wavelength array from photom reference file
        sens_resp : float array
            Response (flux calibration) array from photom reference file
        """
        # Initialize the simulated dispersed image
        self.simulated_image = np.zeros(self.dims, float)

        width = 1.0
        height = 1.0
        x0 = self.xs + 0.5 * width
        y0 = self.ys + 0.5 * height

        # Compute the WCS transforms
        # Setup the transforms we need from the input WCS objects
        sky_to_imgxy = self.grism_wcs.get_transform("world", "detector")
        imgxy_to_grismxy = self.grism_wcs.get_transform("detector", "grism_detector")

        # Get x/y positions in the grism image corresponding to wmin and wmax:
        # Start with RA/Dec of the input pixel position in segmentation map,
        x0_sky, y0_sky = self.seg_wcs(x0, y0)
        # then convert to x/y in the direct image frame corresponding
        # to the grism image,
        x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, 1, order)
        # then finally convert to x/y in the grism image frame
        xwmin, ywmin = imgxy_to_grismxy(x0_xy + self.xoffset, y0_xy + self.yoffset, wmin, order)
        xwmax, ywmax = imgxy_to_grismxy(x0_xy + self.xoffset, y0_xy + self.yoffset, wmax, order)
        dxw = xwmax - xwmin
        dyw = ywmax - ywmin

        # Create list of wavelengths on which to compute dispersed pixels
        lams = np.array([self.pivlam] * len(self.fluxes))
        dw = np.abs((wmax - wmin) / (dyw - dxw))
        dlam = np.median(
            dw / self.oversample_factor
        )  # TODO: validate that just taking median is ok
        lambdas = np.arange(wmin, wmax + dlam, dlam)
        n_lam = len(lambdas)

        # Compute lists of x/y positions in the grism image for
        # the set of desired wavelengths:
        # x/y in image frame of grism image is the same for all wavelengths
        x0_sky = np.repeat(x0_sky[np.newaxis, :], n_lam, axis=0)
        y0_sky = np.repeat(y0_sky[np.newaxis, :], n_lam, axis=0)
        source_ids = np.repeat(self.source_ids_per_pixel[np.newaxis, :], n_lam, axis=0)
        fluxes = np.repeat(self.fluxes[np.newaxis, :], n_lam, axis=0)
        x0_xy, y0_xy, _, _ = sky_to_imgxy(x0_sky, y0_sky, lambdas, order)

        # Convert to x/y in grism frame.
        # lambdas needs same shape as x0_xy to be indexed by np.take below
        lambdas = np.repeat(lambdas[:, np.newaxis], x0_xy.shape[1], axis=1)
        x0s, y0s = imgxy_to_grismxy(x0_xy + self.xoffset, y0_xy + self.yoffset, lambdas, order)
        # x0s, y0s now have shape (n_lam, n_pixels)

        # Compute arrays of dispersed pixel locations and areas
        padding = 1
        naxis = self.dims[::-1]
        # If none of the dispersed pixel indexes are within the image frame,
        # return a null result without wasting time doing other computations
        if x0s.min() >= naxis[0] or x0s.max() < 0 or y0s.min() >= naxis[1] or y0s.max() < 0:
            log.info(f"No dispersed pixels within image frame for order {order}.")
            return

        xs, ys, areas, index = get_clipped_pixels(
            x0s, y0s, padding, naxis[0], naxis[1], width, height
        )
        lams = np.take(lambdas, index)
        fluxes = np.take(fluxes, index)

        # compute 1D sensitivity array corresponding to list of wavelengths
        # TODO: what wavelength unit does this expect?
        sens, no_cal = create_1d_sens(lams, sens_waves, sens_resp)

        # Compute countrates for dispersed pixels. Note that dispersed pixel
        # values are naturally in units of physical fluxes, so we divide out
        # the sensitivity (flux calibration) values to convert to units of
        # countrate (DN/s).
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero")
            counts = fluxes * lams * areas / (sens * self.oversample_factor)
        counts[no_cal] = 0.0  # set to zero where no flux cal info available

        # keep track of source IDs for each dispersed pixel
        dispersed_source_ids = np.take(source_ids, index)

        # Build the dispersed image and make a slit model for each source
        for this_sid in self.source_ids:
            this_sid_idx = dispersed_source_ids == this_sid
            this_xs = xs[this_sid_idx]
            this_ys = ys[this_sid_idx]
            this_flxs = counts[this_sid_idx]
            if len(this_xs) == 0:
                continue

            img, bounds = self._build_dispersed_image_of_source(this_xs, this_ys, this_flxs)
            self.simulated_image[bounds[2] : bounds[3] + 1, bounds[0] : bounds[1] + 1] += img

            slit = self.construct_slitmodel(img, bounds, this_sid, order)
            self.simul_slits.slits.append(slit)
            self.simul_slits_order.append(order)
            self.simul_slits_sid.append(this_sid)

    @staticmethod
    def _build_dispersed_image_of_source(x, y, flux):
        minx = int(min(x))
        maxx = int(max(x))
        miny = int(min(y))
        maxy = int(max(y))
        a = sparse.coo_matrix(
            (flux, (y - miny, x - minx)), shape=(maxy - miny + 1, maxx - minx + 1)
        ).toarray()
        bounds = [minx, maxx, miny, maxy]
        return a, bounds

    @staticmethod
    def construct_slitmodel(
        chunk_data,
        bounds,
        sid,
        order,
    ):
        """
        Turn output image from a chunk into a slit model.

        Parameters
        ----------
        chunk_data : np.ndarray
            Dispersed model of segmentation map source
        bounds : list
            The bounds of the object
        sid : int
            The source ID
        order : int
            The spectral order number

        Returns
        -------
        slit : `jwst.datamodels.SlitModel`
            Slit model containing the dispersed pixel values
        """
        [thisobj_minx, thisobj_maxx, thisobj_miny, thisobj_maxy] = bounds
        slit = datamodels.SlitModel()
        slit.source_id = sid
        slit.name = f"source_{sid}"
        slit.xstart = thisobj_minx
        slit.xsize = thisobj_maxx - thisobj_minx + 1
        slit.ystart = thisobj_miny
        slit.ysize = thisobj_maxy - thisobj_miny + 1
        slit.meta.wcsinfo.spectral_order = order
        slit.data = chunk_data[thisobj_miny : thisobj_maxy + 1, thisobj_minx : thisobj_maxx + 1]

        return slit
