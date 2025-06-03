import time
import numpy as np
import multiprocessing as mp

from scipy import sparse

from stdatamodels.jwst import datamodels

from .disperse import dispersed_pixel

import logging

from photutils.background import Background2D, MedianBackground
from astropy.stats import SigmaClip

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def _disperse_multiprocess(pars, max_cpu):
    ctx = mp.get_context("spawn")
    with ctx.Pool(max_cpu) as mypool:
        all_res = mypool.starmap(dispersed_pixel, pars)

    return all_res


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
        direct_images,
        segmap_model,
        grism_wcs,
        filter_name,
        source_id=None,
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
        self.dir_image_names = direct_images
        self.seg = segmap_model.data
        all_ids = np.array(list(set(np.ravel(self.seg))))
        self.source_ids = _select_ids(source_id, all_ids)
        self.filter = filter_name
        self.sed_file = sed_file  # should always be NONE for baseline pipeline (use flat SED)
        self.renormalize = renormalize
        self.max_cpu = max_cpu
        self.xoffset = offsets[0]
        self.yoffset = offsets[1]

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
        self.create_pixel_list()

        # Initialize the list of slits
        self.simul_slits = datamodels.MultiSlitModel()
        self.simul_slits_order = []
        self.simul_slits_sid = []

    def create_pixel_list(self):
        """
        Create a list of pixels to be dispersed, grouped per object ID.

        When ID is None, all sources in the segmentation map are processed.
        """
        self.xs = []
        self.ys = []
        for source_id in self.source_ids:
            ys, xs = np.nonzero(self.seg == source_id)
            if len(xs) > 0 and len(ys) > 0:
                self.xs.append(xs)
                self.ys.append(ys)

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

    def disperse_all(
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

        # Loop over all source ID's from segmentation map
        pool_args = []
        for i in range(len(self.source_ids)):
            disperse_chunk_args = [
                i,
                order,
                wmin,
                wmax,
                sens_waves,
                sens_resp,
            ]
            pool_args.append(disperse_chunk_args)

        t0 = time.time()
        if self.max_cpu > 1:
            # put this log message here to avoid printing it for every chunk
            log.info(f"Using multiprocessing with {self.max_cpu} cores to compute dispersion")

        disperse_chunk_output = []
        for i in range(len(self.source_ids)):
            disperse_chunk_output.append(self.disperse_chunk(*pool_args[i]))
        t1 = time.time()
        log.info(f"Wall clock time for disperse_chunk order {order}: {(t1 - t0):.1f} sec")

        # Collect results into simulated image and slit models
        for this_output in disperse_chunk_output:
            [this_image, this_bounds, this_sid, this_order] = this_output
            slit = self.construct_slitmodel_for_chunk(this_image, this_bounds, this_sid, this_order)
            self.simulated_image += this_image
            if slit is not None:
                self.simul_slits.slits.append(slit)
                self.simul_slits_order.append(this_order)
                self.simul_slits_sid.append(this_sid)

    def disperse_chunk(
        self,
        c,
        order,
        wmin,
        wmax,
        sens_waves,
        sens_resp,
    ):
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
        this_object : np.ndarray
            2D array of dispersed pixel values for the source
        thisobj_bounds : list
            [minx, maxx, miny, maxy] bounds of the object
        sid : int
            Source ID
        order : int
            Spectral order number
        """
        sid = int(self.source_ids[c])
        self.order = order
        self.wmin = wmin
        self.wmax = wmax
        self.sens_waves = sens_waves
        self.sens_resp = sens_resp
        log.info(
            f"Dispersing source {sid}, order {self.order}. "
            f"Source contains {len(self.xs[c])} pixels."
        )

        # Loop over all pixels in list for object "c"
        pars = []  # initialize params for this object
        for i in range(len(self.xs[c])):
            # Here "i" is just an index into the pixel list for the object
            # being processed, as opposed to the ID number of the object itself
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
            pars_i = [
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
                i,
                self.dims[::-1],
                2,
                self.extrapolate_sed,
                self.xoffset,
                self.yoffset,
            ]
            pars.append(pars_i)
            # if i == 0:
            #    print([type(arg) for arg in pars_i]) #all these need to be pickle-able

        # pass parameters into dispersed_pixel, either using multiprocessing or not
        time1 = time.time()
        if self.max_cpu > 1:
            all_res = _disperse_multiprocess(pars, self.max_cpu)
        else:
            all_res = []
            for i in range(len(pars)):
                all_res.append(dispersed_pixel(*pars[i]))

        # Initialize blank image for this source
        this_object = np.zeros(self.dims, float)
        nres = 0
        bounds = []
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
            this_object[miny : maxy + 1, minx : maxx + 1] += a
            bounds.append([minx, maxx, miny, maxy])

        time2 = time.time()
        log.debug(f"Elapsed time {time2 - time1} sec")
        # figure out global bounds of object
        if len(bounds) > 0:
            bounds = np.array(bounds)
            thisobj_minx = int(np.min(bounds[:, 0]))
            thisobj_maxx = int(np.max(bounds[:, 1]))
            thisobj_miny = int(np.min(bounds[:, 2]))
            thisobj_maxy = int(np.max(bounds[:, 3]))
            thisobj_bounds = [thisobj_minx, thisobj_maxx, thisobj_miny, thisobj_maxy]
            return (this_object, thisobj_bounds, sid, order)
        return (this_object, None, sid, order)

    @staticmethod
    def construct_slitmodel_for_chunk(
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
        if bounds is None:
            return None

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
