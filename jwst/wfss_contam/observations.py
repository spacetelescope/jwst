import logging
import multiprocessing as mp
import time

import numpy as np
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from stdatamodels.jwst import datamodels

from jwst.wfss_contam.disperse import disperse

log = logging.getLogger(__name__)

__all__ = ["background_subtract", "Observation"]


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
    data : ndarray
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
    data : ndarray
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
    all_ids : ndarray
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
    """
    Define an observation leading to a single grism image.

    The Observation class is responsible for calling the various WCS transforms that convert
    a direct image and a segmentation image into a simulation of the grism image, making
    assumptions about the spectral properties of the direct image sources.
    When the disperse_order method is called one or more times, two products are created:
    the simulated dispersed image (simulated_image attribute) and
    the simulated MultiSlitModel (simulated_slits attribute).
    """

    def __init__(
        self,
        direct_image,
        segmentation_map,
        grism_wcs,
        direct_image_wcs,
        boundaries=None,
        max_cpu=1,
        max_pixels_per_chunk=5e4,
        oversample_factor=2,
        phot_per_lam=True,
    ):
        """
        Initialize all data and metadata for a given observation.

        Parameters
        ----------
        direct_image : np.ndarray
            Direct imaging data.
        segmentation_map : np.ndarray
            Segmentation map data.
        grism_wcs : `~gwcs.wcs.WCS`
            WCS object from grism image
        direct_image_wcs : `~gwcs.wcs.WCS`
            WCS object from direct image
        boundaries : list, optional
            Start/Stop coordinates of the FOV within the larger seed image.
        max_cpu : int, optional
            Max number of cpu's to use when multiprocessing
        max_pixels_per_chunk : int, optional
            Maximum number of pixels per chunk when dispersing sources
        oversample_factor : int, optional
            Factor by which to oversample the wavelength grid
        phot_per_lam : bool, optional
            Whether to compute photometry per wavelength bin (True) or per pixel (False).
            This depends on how the photom reference file has been delivered.
            True should be used for NIRCam, and False should be used for NIRISS.
        """
        if boundaries is None:
            boundaries = []
        # Load all the info for this grism mode
        self.direct_image_wcs = direct_image_wcs
        self.grism_wcs = grism_wcs
        self.seg = segmentation_map
        all_ids = list(set(np.ravel(self.seg)))
        all_ids.remove(0)  # Remove the background ID
        self.source_ids = all_ids
        self.max_cpu = max_cpu
        self.max_pixels_per_chunk = max_pixels_per_chunk
        self.oversample_factor = oversample_factor
        self.phot_per_lam = phot_per_lam

        # ensure the direct image has background subtracted
        self.dimage = background_subtract(direct_image)

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
        self.naxis = self.dims[::-1]

        # Create lists of pixels labeled in segmentation map
        self._create_pixel_list()

        # Initialize the output MultiSlitModel
        self.simulated_slits = datamodels.MultiSlitModel()

        # Initialize the simulated dispersed image
        self.simulated_image = np.zeros(self.dims, float)

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

    def chunk_sources(
        self, order, wmin, wmax, sens_waves, sens_response, selected_ids=None, max_pixels=1e5
    ):
        """
        Chunk the sources into groups of max_pixels.

        Parameters
        ----------
        max_pixels : int, optional
            Maximum number of pixels per chunk.

        Returns
        -------
        disperse_args : list[list]
            Outer list has length number of groups, and each inner list contains
            the arguments to disperse() for that group
            in the format that multiprocessing starmap expects.
        """
        chunks = []
        current_chunk = []
        current_size = 0

        source_ids = _select_ids(selected_ids, self.source_ids)

        for sid in source_ids:
            n_pixels = np.sum(self.seg == sid)
            if n_pixels > max_pixels:
                log.warning(
                    f"Source {sid} has {n_pixels} pixels, which exceeds the maximum number "
                    f"of pixels per chunk ({max_pixels}). Skipping this source. "
                    "Consider increasing max_pixels, and/or check to ensure the segmentation map "
                    "looks reasonable for that source."
                )
                continue
            if current_size + n_pixels > max_pixels:
                chunks.append(current_chunk)
                current_chunk = []
                current_size = 0
            current_chunk.append(sid)
            current_size += n_pixels

        if current_chunk:
            chunks.append(current_chunk)

        disperse_args = []
        for source_ids in chunks:
            isin = np.isin(self.source_ids_per_pixel, source_ids)
            source_ids_per_pixel = self.source_ids_per_pixel[isin]
            xs = self.xs[isin]
            ys = self.ys[isin]
            fluxes = self.fluxes[isin]
            disperse_args.append(
                [
                    xs,
                    ys,
                    fluxes,
                    source_ids_per_pixel,
                    order,
                    wmin,
                    wmax,
                    sens_waves,
                    sens_response,
                    self.direct_image_wcs,
                    self.grism_wcs,
                    self.naxis,
                    self.oversample_factor,
                    self.phot_per_lam,
                ]
            )

        return disperse_args

    def disperse_order(self, order, wmin, wmax, sens_waves, sens_response, selected_ids=None):
        """
        Disperse the sources for a given spectral order, with multiprocessing.

        The simulated_slits and simulated_image attributes are updated in place.

        Parameters
        ----------
        order : int
            Spectral order to process
        wmin : float
            Minimum wavelength for dispersed spectra
        wmax : float
            Maximum wavelength for dispersed spectra
        sens_waves : ndarray
            Wavelength array from photom reference file
        sens_response : ndarray
            Response (flux calibration) array from photom reference file
        selected_ids : list, optional
            List of source IDs to process. If None, all sources are processed.
        """
        # generate lists of input parameters for the disperse function
        # for each chunk of sources
        disperse_args = self.chunk_sources(
            order,
            wmin,
            wmax,
            sens_waves,
            sens_response,
            selected_ids=selected_ids,
            max_pixels=self.max_pixels_per_chunk,
        )
        t0 = time.time()
        if self.max_cpu > 1:
            # Use multiprocessing to disperse the sources
            log.info(
                f"Using {self.max_cpu} CPU cores for multiprocessing "
                f"{len(self.source_ids)} sources in {len(disperse_args)} chunks."
            )
            ctx = mp.get_context("spawn")
            with ctx.Pool(self.max_cpu) as mypool:
                all_res = mypool.starmap(disperse, disperse_args)
        else:
            all_res = [disperse(*args) for args in disperse_args]
        t1 = time.time()
        log.info(f"Wall clock time for disperse_chunk order {order}: {(t1 - t0):.1f} sec")

        # Combine the results from all chunks
        for results in all_res:
            if results is None:
                # None of the sources in this chunk for this order had pixels on the detector
                continue
            for sid in results:
                bounds = results[sid]["bounds"]
                img = results[sid]["image"]
                slit = _construct_slitmodel(img, bounds, sid, order)
                self.simulated_image[bounds[2] : bounds[3] + 1, bounds[0] : bounds[1] + 1] += img
                self.simulated_slits.slits.append(slit)


def _construct_slitmodel(
    img,
    bounds,
    sid,
    order,
):
    """
    Turn an output image from a single source/order into a SlitModel.

    Parameters
    ----------
    img : ndarray
        Dispersed model image of segmentation map source
    bounds : list
        The bounds of the object in relation to the full-frame image.
    sid : int
        The source ID
    order : int
        The spectral order

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
    slit.data = img

    return slit
