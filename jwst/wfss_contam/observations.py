import logging
import multiprocessing as mp
import time
import warnings

import numpy as np
from astropy.stats import SigmaClip
from astropy.utils.exceptions import AstropyUserWarning
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
    with warnings.catch_warnings():
        # there can be multiple different AstropyUserWarning messages here about NaN and Inf values
        warnings.filterwarnings("ignore", category=AstropyUserWarning)
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
        source_ids = _select_ids(selected_ids, self.source_ids)
        max_pixels = int(max_pixels)

        # Create a mask for selected sources
        selected_mask = np.isin(self.source_ids_per_pixel, source_ids)

        # Get pixels for selected sources
        selected_xs = self.xs[selected_mask]
        selected_ys = self.ys[selected_mask]
        selected_fluxes = self.fluxes[selected_mask]
        selected_source_ids = self.source_ids_per_pixel[selected_mask]

        # Sort by source ID to keep sources mostly together
        # This reduces the number of times we have to call build_dispersed_image_of_source
        # within disperse()
        sort_indices = np.argsort(selected_source_ids)
        sorted_xs = selected_xs[sort_indices]
        sorted_ys = selected_ys[sort_indices]
        sorted_fluxes = selected_fluxes[sort_indices]
        sorted_source_ids = selected_source_ids[sort_indices]

        # Split into chunks of max_pixels
        total_pixels = len(sorted_xs)
        n_chunks = int(np.ceil(total_pixels / max_pixels))

        log.info(
            f"Splitting {total_pixels} pixels from {len(source_ids)} sources into {n_chunks} chunks"
        )

        disperse_args = []
        for i in range(n_chunks):
            start_idx = i * max_pixels
            end_idx = min((i + 1) * max_pixels, total_pixels)

            chunk_xs = sorted_xs[start_idx:end_idx]
            chunk_ys = sorted_ys[start_idx:end_idx]
            chunk_fluxes = sorted_fluxes[start_idx:end_idx]
            chunk_source_ids = sorted_source_ids[start_idx:end_idx]

            disperse_args.append(
                [
                    chunk_xs,
                    chunk_ys,
                    chunk_fluxes,
                    chunk_source_ids,
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
            pool = ctx.Pool(self.max_cpu)
            try:
                all_res = pool.starmap(disperse, disperse_args)
            except Exception as e:
                log.error(f"Error during parallel processing: {e}")
                raise
            finally:
                pool.close()
                pool.join()
        else:
            all_res = [disperse(*args) for args in disperse_args]
        t1 = time.time()
        log.info(f"Wall clock time for disperse_chunk order {order}: {(t1 - t0):.1f} sec")

        # Combine results from all chunks, aggregating by source ID
        source_results = {}
        for results in all_res:
            if results is None:
                # None of the sources in this chunk for this order had pixels on the detector
                continue
            for sid in results:
                _aggregate_by_source(results, sid, source_results)

        # Now add the combined results to the simulation
        for sid in source_results:
            bounds = source_results[sid]["bounds"]
            img = source_results[sid]["image"]
            slit = _construct_slitmodel(img, bounds, sid, order)
            self.simulated_image[bounds[2] : bounds[3] + 1, bounds[0] : bounds[1] + 1] += img
            self.simulated_slits.slits.append(slit)


def _aggregate_by_source(results, sid, source_results):
    if sid not in source_results:
        source_results[sid] = {
            "bounds": results[sid]["bounds"],
            "image": results[sid]["image"].copy(),
        }
    else:
        # Combine bounds
        old_bounds = source_results[sid]["bounds"]
        new_bounds = results[sid]["bounds"]
        combined_bounds = [
            min(old_bounds[0], new_bounds[0]),
            max(old_bounds[1], new_bounds[1]),
            min(old_bounds[2], new_bounds[2]),
            max(old_bounds[3], new_bounds[3]),
        ]

        # Create combined image with the union of bounds
        combined_shape = (
            combined_bounds[3] - combined_bounds[2] + 1,
            combined_bounds[1] - combined_bounds[0] + 1,
        )
        combined_image = np.zeros(combined_shape, dtype=float)

        # Add existing image to combined image
        old_y_start = old_bounds[2] - combined_bounds[2]
        old_y_end = old_y_start + source_results[sid]["image"].shape[0]
        old_x_start = old_bounds[0] - combined_bounds[0]
        old_x_end = old_x_start + source_results[sid]["image"].shape[1]
        combined_image[old_y_start:old_y_end, old_x_start:old_x_end] += source_results[sid]["image"]

        # Add new image to combined image
        new_y_start = new_bounds[2] - combined_bounds[2]
        new_y_end = new_y_start + results[sid]["image"].shape[0]
        new_x_start = new_bounds[0] - combined_bounds[0]
        new_x_end = new_x_start + results[sid]["image"].shape[1]
        combined_image[new_y_start:new_y_end, new_x_start:new_x_end] += results[sid]["image"]

        # Update source results
        source_results[sid] = {"bounds": combined_bounds, "image": combined_image}


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
