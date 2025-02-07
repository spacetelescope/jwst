import time
import multiprocessing
import numpy as np

from scipy import sparse

from stdatamodels.jwst import datamodels

from .disperse import dispersed_pixel

import logging

from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def background_subtract(
    data, box_size=None, filter_size=(3, 3), sigma=3.0, exclude_percentile=30.0
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


class Observation:
    """Define an observation leading to a single grism image."""

    def __init__(
        self,
        direct_images,
        segmap_model,
        grism_wcs,
        filter_name,
        source_id=0,
        sed_file=None,
        extrapolate_sed=False,
        boundaries=None,
        offsets=None,
        renormalize=True,
        max_cpu=1,
    ):
        """
        Initialize all data and metadata for a given observation.

        Creates lists of
        direct image pixel values for selected objects.

        Parameters
        ----------
        direct_images : List of strings
            List of file name(s) containing direct imaging data
        segmap_model : `jwst.datamodels.ImageModel`
            Segmentation map model
        grism_wcs : gwcs object
            WCS object from grism image
        filter_name : str
            Filter name
        source_id : int, optional, default 0
            ID of source to process. If 0, all sources processed.
        sed_file : str, optional, default None
            Name of Spectral Energy Distribution (SED) file containing datasets matching
            the ID in the segmentation file and each consisting of a [[lambda],[flux]] array.
        extrapolate_sed : bool, optional, default False
            Flag indicating whether to extrapolate wavelength range of SED
        boundaries : list, optional, default []
            Start/Stop coordinates of the FOV within the larger seed image.
        offsets : list, optional, default [0,0]
            Offset values for x and y axes
        renormalize : bool, optional, default True
            Flag indicating whether to renormalize SED's
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
        self.source_id = source_id
        self.source_ids = []
        self.dir_image_names = direct_images
        self.seg = segmap_model.data
        self.filter = filter_name
        self.sed_file = sed_file  # should always be NONE for baseline pipeline (use flat SED)
        self.cache = False
        self.renormalize = renormalize
        self.max_cpu = max_cpu
        self.xoffset = offsets[0]
        self.yoffset = offsets[1]

        # Set the limits of the dispersed image to be simulated
        if len(boundaries) == 0:
            log.debug("No boundaries passed.")
            self.xstart = 0
            self.xend = self.xstart + self.dims[0] - 1
            self.ystart = 0
            self.yend = self.ystart + self.dims[1] - 1
        else:
            self.xstart, self.xend, self.ystart, self.yend = boundaries
        self.dims = (self.yend - self.ystart + 1, self.xend - self.xstart + 1)
        log.debug(f"Using simulated image size of {self.dims[1]} {self.dims[0]}")

        # Allow for SED extrapolation
        self.extrapolate_sed = extrapolate_sed
        if self.extrapolate_sed:
            log.warning("SED Extrapolation turned on.")

        # Create pixel lists for sources labeled in segmentation map
        self.create_pixel_list()

    def create_pixel_list(self):
        """Create a list of pixels to be dispersed, grouped per object ID."""
        if self.source_id == 0:
            # When source_id=0, all sources in the segmentation map are processed.
            # This creates a huge list of all x,y pixel indices that have non-zero values
            # in the seg map, sorted by those indices belonging to a particular source ID.
            self.xs = []
            self.ys = []
            all_ids = np.array(list(set(np.ravel(self.seg))))
            all_ids = all_ids[all_ids > 0]
            self.source_ids = all_ids
            log.info(f"Loading {len(all_ids)} sources from segmentation map")
            for source_id in all_ids:
                ys, xs = np.nonzero(self.seg == source_id)
                if len(xs) > 0 and len(ys) > 0:
                    self.xs.append(xs)
                    self.ys.append(ys)

        else:
            # Process only the given source ID
            log.info(f"Loading source {self.source_id} from segmentation map")
            ys, xs = np.nonzero(self.seg == self.source_id)
            if len(xs) > 0 and len(ys) > 0:
                self.xs = [xs]
                self.ys = [ys]
                self.source_ids = [self.source_id]

        # Populate lists of direct image flux values for the sources.
        self.fluxes = {}
        for dir_image_name in self.dir_image_names:
            log.info(f"Using direct image {dir_image_name}")
            with datamodels.open(dir_image_name) as model:
                dimage = model.data
                dimage = background_subtract(dimage)

                if self.sed_file is None:
                    # Default pipeline will use sed_file=None, so we need to compute
                    # photometry values that used to come from HST-style header keywords.
                    # Set pivlam, in units of microns, based on filter name.
                    pivlam = float(self.filter[1:4]) / 100.0

                    # Use pixel fluxes from the direct image.
                    self.fluxes[pivlam] = []
                    for i in range(len(self.source_ids)):
                        # This loads lists of pixel flux values for each source
                        # from the direct image
                        self.fluxes[pivlam].append(dimage[self.ys[i], self.xs[i]])

                else:
                    # Use an SED file. Need to normalize the object stamps.
                    for source_id in self.source_ids:
                        vg = self.seg == source_id
                        dnew = dimage
                        if self.renormalize:
                            sum_seg = np.sum(dimage[vg])  # But normalize by the whole flux
                            if sum_seg != 0:
                                dimage[vg] /= sum_seg
                        else:
                            log.debug("not renormalizing sources to unity")

                    self.fluxes["sed"] = []
                    for i in range(len(self.source_ids)):
                        self.fluxes["sed"].append(dnew[self.ys[i], self.xs[i]])

    def disperse_all(self, order, wmin, wmax, sens_waves, sens_resp, cache=False):
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
        if cache:
            log.debug("Object caching ON")
            self.cache = True
            self.cached_object = {}

        # Initialize the simulated dispersed image
        self.simulated_image = np.zeros(self.dims, float)

        # Loop over all source IDs from segmentation map
        for i in range(len(self.source_ids)):
            if self.cache:
                self.cached_object[i] = {}
                self.cached_object[i]["x"] = []
                self.cached_object[i]["y"] = []
                self.cached_object[i]["f"] = []
                self.cached_object[i]["w"] = []
                self.cached_object[i]["minx"] = []
                self.cached_object[i]["maxx"] = []
                self.cached_object[i]["miny"] = []
                self.cached_object[i]["maxy"] = []

            self.disperse_chunk(i, order, wmin, wmax, sens_waves, sens_resp)

    def disperse_chunk(self, c, order, wmin, wmax, sens_waves, sens_resp):
        """
        Compute dispersion for a single source; to be called after create_pixel_list().

        Parameters
        ----------
        c : int
            Chunk (source) number to process
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

        Returns
        -------
        np.ndarray
            2D dispersed image for this source
        """
        sid = int(self.source_ids[c])
        self.order = order
        self.wmin = wmin
        self.wmax = wmax
        self.sens_waves = sens_waves
        self.sens_resp = sens_resp
        log.info(f"Dispersing source {sid}, order {self.order}")
        pars = []  # initialize params for this object

        # Loop over all pixels in list for object "c"
        log.debug(f"source contains {len(self.xs[c])} pixels")
        for i in range(len(self.xs[c])):
            # Here "i" just indexes the pixel list for the object being processed

            # xc, yc are the coordinates of the central pixel of the group
            # of pixels surrounding the direct image pixel index
            width = 1.0
            height = 1.0
            xc = self.xs[c][i] + 0.5 * width
            yc = self.ys[c][i] + 0.5 * height

            # "lams" is the array of wavelengths previously stored in flux list
            # and correspond to the central wavelengths of the filters used in
            # the input direct image(s). For the simple case of 1 combined direct image,
            # this contains a single value (e.g. 4.44 for F444W).

            # "fluxes" is the array of pixel values from the direct image(s).
            # For the simple case of 1 combined direct image, this contains a
            # a single value (just like "lams").
            fluxes, lams = map(
                np.array,
                zip(
                    *[
                        (self.fluxes[lm][c][i], lm)
                        for lm in sorted(self.fluxes.keys())
                        if self.fluxes[lm][c][i] != 0
                    ],
                    strict=True,
                ),
            )

            pars_i = (
                xc,
                yc,
                width,
                height,
                lams,
                fluxes,
                self.order,
                self.wmin,
                self.wmax,
                self.sens_waves,
                self.sens_resp,
                self.seg_wcs,
                self.grism_wcs,
                i,  # TODO: this is not the source_id as the docstring to dispersed_pixel says
                self.dims[::-1],
                2,
                self.extrapolate_sed,
                self.xoffset,
                self.yoffset,
            )

            pars.append(pars_i)
            # now have full pars list for all pixels for this object

        time1 = time.time()
        if self.max_cpu > 1:
            ctx = multiprocessing.get_context("forkserver")
            mypool = ctx.Pool(self.max_cpu)  # Create the pool
            all_res = mypool.imap_unordered(dispersed_pixel, pars)  # Fill the pool
            mypool.close()  # Drain the pool
        else:
            all_res = []
            for i in range(len(pars)):
                all_res.append(dispersed_pixel(*pars[i]))

        # Initialize blank image for this source
        this_object = np.zeros(self.dims, float)

        nres = 0
        for pp in all_res:
            if pp is None:
                continue

            nres += 1
            x, y, _, w, f, *_ = pp

            # skip results that don't have pixels in the field
            if len(x) < 1:
                continue

            minx = int(min(x))
            maxx = int(max(x))
            miny = int(min(y))
            maxy = int(max(y))
            a = sparse.coo_matrix(
                (f, (y - miny, x - minx)), shape=(maxy - miny + 1, maxx - minx + 1)
            ).toarray()

            # Accumulate results into simulated images
            self.simulated_image[miny : maxy + 1, minx : maxx + 1] += a
            this_object[miny : maxy + 1, minx : maxx + 1] += a

            if self.cache:
                self.cached_object[c]["x"].append(x)
                self.cached_object[c]["y"].append(y)
                self.cached_object[c]["f"].append(f)
                self.cached_object[c]["w"].append(w)
                self.cached_object[c]["minx"].append(minx)
                self.cached_object[c]["maxx"].append(maxx)
                self.cached_object[c]["miny"].append(miny)
                self.cached_object[c]["maxy"].append(maxy)

        time2 = time.time()
        log.debug(f"Elapsed time {time2 - time1} sec")

        return this_object

    def disperse_all_from_cache(self, trans=None):
        """
        Compute dispersed pixel values for all sources identified in the segmentation map.

        Load data from cache where available. Currently not used.

        Parameters
        ----------
        trans : function
            Transmission function to apply to the flux values

        Returns
        -------
        np.ndarray
            2D dispersed image for this source

        Notes
        -----
        The return value of `this_object` appears to be a bug.
        However, this is currently not used, and if the INS team wants to re-enable
        caching, all functions here need updating anyway, so not fixing at this time.
        """
        if not self.cache:
            return

        self.simulated_image = np.zeros(self.dims, float)

        for i in range(len(self.source_ids)):
            this_object = self.disperse_chunk_from_cache(i, trans=trans)

        return this_object

    def disperse_chunk_from_cache(self, c, trans=None):
        """
        Compute dispersion for a single source; to be called after create_pixel_list().

        Load data from cache where available. Currently not used.

        Parameters
        ----------
        c : int
            Chunk (source) number to process
        trans : function
            Transmission function to apply to the flux values

        Returns
        -------
        np.ndarray
            2D dispersed image for this source
        """
        if not self.cache:
            return

        time1 = time.time()

        # Initialize blank image for this object
        this_object = np.zeros(self.dims, float)

        if trans is not None:
            log.debug("Applying a transmission function...")

        for i in range(len(self.cached_object[c]["x"])):
            x = self.cached_object[c]["x"][i]
            y = self.cached_object[c]["y"][i]
            f = self.cached_object[c]["f"][i] * 1.0
            w = self.cached_object[c]["w"][i]

            if trans is not None:
                f *= trans(w)

            minx = self.cached_object[c]["minx"][i]
            maxx = self.cached_object[c]["maxx"][i]
            miny = self.cached_object[c]["miny"][i]
            maxy = self.cached_object[c]["maxy"][i]

            a = sparse.coo_matrix(
                (f, (y - miny, x - minx)), shape=(maxy - miny + 1, maxx - minx + 1)
            ).toarray()

            # Accumulate the results into the simulated images
            self.simulated_image[miny : maxy + 1, minx : maxx + 1] += a
            this_object[miny : maxy + 1, minx : maxx + 1] += a

        time2 = time.time()
        log.debug(f"Elapsed time {time2 - time1} sec")

        return this_object
