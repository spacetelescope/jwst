import time
import numpy as np
from multiprocessing import Pool

from scipy import sparse

from stdatamodels.jwst import datamodels

from .disperse import dispersed_pixel

import logging
import warnings

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class Observation:
    """This class defines an actual observation. It is tied to a single grism image."""

    def __init__(self, direct_images, segmap_model, grism_wcs, filter, ID=0,
                 sed_file=None, extrapolate_sed=False,
                 boundaries=[], offsets=[0, 0], renormalize=True, max_cpu=1):

        """
        Initialize all data and metadata for a given observation. Creates lists of
        direct image pixel values for selected objects.

        Parameters
        ----------
        direct_images : List of strings
            List of file name(s) containing direct imaging data
        segmap_model : `jwst.datamodels.ImageModel`
            Segmentation map model
        grism_wcs : gwcs object
            WCS object from grism image
        filter : str
            Filter name
        ID : int or list-like, optional
            ID(s) of source to process. If zero, all sources processed.
        sed_file : str
            Name of Spectral Energy Distribution (SED) file containing datasets matching
            the ID in the segmentation file and each consisting of a [[lambda],[flux]] array.
        extrapolate_sed : bool
            Flag indicating whether to extrapolate wavelength range of SED
        boundaries : tuple
            Start/Stop coordinates of the FOV within the larger seed image.
        renormalize : bool
            Flag indicating whether to renormalize SED's
        max_cpu : int
            Max number of cpu's to use when multiprocessing
        """

        # Load all the info for this grism mode
        self.seg_wcs = segmap_model.meta.wcs
        self.grism_wcs = grism_wcs
        self.ID = ID
        self.IDs = []
        self.dir_image_names = direct_images
        self.seg = segmap_model.data
        self.filter = filter
        self.sed_file = sed_file   # should always be NONE for baseline pipeline (use flat SED)
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

        # Initialize the list of slits
        self.simul_slits = datamodels.MultiSlitModel()
        self.simul_slits_order = []
        self.simul_slits_sid = []

    def create_pixel_list(self):
        # Create a list of pixels to be dispersed, grouped per object ID.

        if self.ID == 0:
            # When ID=0, all sources in the segmentation map are processed.
            # This creates a huge list of all x,y pixel indices that have non-zero values
            # in the seg map, sorted by those indices belonging to a particular source ID.
            all_IDs = np.array(list(set(np.ravel(self.seg))))
            all_IDs = all_IDs[all_IDs > 0]
            self.IDs = all_IDs
            log.info(f"Loading {len(all_IDs)} sources from segmentation map")
        elif isinstance(self.ID, int):
            # Process only the given source ID
            log.info(f"Loading source {self.ID} from segmentation map")
            self.IDs = [self.ID]
        elif isinstance(self.ID, (list, np.array)):
            # Process only the given list of source IDs
            log.info(f"Loading {len(self.ID)} of {len(list(set(np.ravel(self.seg))))} selected sources from segmentation map")
            self.IDs = self.ID
        else:
            raise ValueError("ID must be an integer or a list of integers")

        self.xs = []
        self.ys = []
        for ID in self.IDs:
            ys, xs = np.nonzero(self.seg == ID)
            if len(xs) > 0 and len(ys) > 0:
                self.xs.append(xs)
                self.ys.append(ys)

        print("length of xs and ys", len(self.xs), len(self.ys))

        # Populate lists of direct image flux values for the sources.
        self.fluxes = {}
        for dir_image_name in self.dir_image_names:

            log.info(f"Using direct image {dir_image_name}")
            with datamodels.open(dir_image_name) as model:
                dimage = model.data

                if self.sed_file is None:
                    # Default pipeline will use sed_file=None, so we need to compute
                    # photometry values that used to come from HST-style header keywords.
                    # Set pivlam, in units of microns, based on filter name.
                    pivlam = float(self.filter[1:4]) / 100.

                    # Use pixel fluxes from the direct image.
                    self.fluxes[pivlam] = []
                    for i in range(len(self.IDs)):
                        # This loads lists of pixel flux values for each source
                        # from the direct image
                        self.fluxes[pivlam].append(dimage[self.ys[i], self.xs[i]])

                else:
                    # Use an SED file. Need to normalize the object stamps.
                    for ID in self.IDs:
                        vg = self.seg == ID
                        dnew = dimage
                        if self.renormalize:
                            sum_seg = np.sum(dimage[vg])  # But normalize by the whole flux
                            if sum_seg != 0:
                                dimage[vg] /= sum_seg
                        else:
                            log.debug("not renormalizing sources to unity")

                    self.fluxes["sed"] = []
                    for i in range(len(self.IDs)):
                        self.fluxes["sed"].append(dnew[self.ys[i], self.xs[i]])

    def disperse_all(self, order, wmin, wmax, cache=False):
        """
        Compute dispersed pixel values for all sources identified in
        the segmentation map.

        Parameters
        ----------
        order : int
            Spectral order number to process
        wmin : float
            Minimum wavelength for dispersed spectra
        wmax : float
            Maximum wavelength for dispersed spectra
        """
        if cache:
            log.debug("Object caching ON")
            self.cache = True
            self.cached_object = {}

        # Initialize the simulated dispersed image
        self.simulated_image = np.zeros(self.dims, float)

        # Loop over all source ID's from segmentation map
        pool_args = []
        for i in range(len(self.IDs)):
            if self.cache:
                self.cached_object[i] = {}
                self.cached_object[i]['x'] = []
                self.cached_object[i]['y'] = []
                self.cached_object[i]['f'] = []
                self.cached_object[i]['w'] = []
                self.cached_object[i]['minx'] = []
                self.cached_object[i]['maxx'] = []
                self.cached_object[i]['miny'] = []
                self.cached_object[i]['maxy'] = []
            disperse_chunk_args = [i, order, wmin, wmax, 
                                   self.IDs[i], self.xs[i], self.ys[i], 
                                   self.fluxes, #check shape!
                                   self.seg_wcs, self.grism_wcs, self.dims, 
                                   self.extrapolate_sed, self.xoffset, self.yoffset]
            pool_args.append(disperse_chunk_args)

        # call disperse_chunk with optional multiprocessing
        t0 = time.time()
        if self.max_cpu > 1:
            log.info(f"Using multiprocessing with {self.max_cpu} cores to compute dispersion")
            with Pool(self.max_cpu) as mypool:
                disperse_chunk_output = mypool.starmap(self.disperse_chunk, pool_args)
        else:
            for i in range(len(self.IDs)):
                disperse_chunk_output = self.disperse_chunk(*pool_args[i])
        t1 = time.time()
        log.info(f"Wall clock time for disperse_chunk order {order}: {(t1-t0):.1f} sec")
        
        # Collect results into simulated image and slit models
        for i, this_output in enumerate(disperse_chunk_output):
            [this_image, this_bounds, this_sid, this_order] = this_output
            slit = self.construct_slitmodel_for_chunk(this_image, this_bounds, this_sid, this_order)
            self.simulated_image += this_image
            if slit is not None:
                self.simul_slits.slits.append(slit)
                self.simul_slits_order.append(this_order)
                self.simul_slits_sid.append(this_sid)

    @staticmethod
    def disperse_chunk(c, order, wmin, wmax, sid, xs, ys, fluxes_dict, seg_wcs, grism_wcs, dims, extrapolate_sed, xoffset, yoffset):
        """
        Method that computes dispersion for a single source.
        To be called after create_pixel_list().
        Static method to enable parallelization

        Parameters
        ----------
        c : int
            Chunk (source) number to process. used to index the fluxes dict
        order : int
            Spectral order to process
        wmin : float
            Minimum wavelength for dispersed spectra
        wmax : float
            Maximum wavelength for dispersed spectra
        sid : int
            Source ID
        xs : np.ndarray
            X-coordinates of the the central pixel of the group of pixels 
            surrounding the direct image pixel index
        ys : np.ndarray
            Y-coordinates of the the central pixel of the group of pixels
            surrounding the direct image pixel index
        fluxes_dict : dict
            Dictionary of fluxes for each direct image.
            fluxes_dict{"lams"} is the array of wavelengths previously stored in flux list
            and correspond to the central wavelengths of the filters used in
            the input direct image(s). For the simple case of 1 combined direct image,
            this contains a single value (e.g. 4.44 for F444W).
            fluxes_dict{"fluxes"} is the array of pixel values from the direct image(s).
            For the simple case of 1 combined direct image, this contains a
            a single value (just like "lams").
        seg_wcs : gwcs object
            WCS object from segmentation map
        grism_wcs : gwcs object
            WCS object from grism image
        dims : tuple
            Dimensions of the grism image
        extrapolate_sed : bool
            Flag indicating whether to extrapolate wavelength range of SED
        xoffset : int
            Pixel offset to apply when computing the dispersion (accounts for offset from source cutout to full frame)
        yoffset : int
            Pixel offset to apply when computing the dispersion (accounts for offset from source cutout to full frame)

        Returns
        -------
        this_object : np.ndarray
            Dispersed model of segmentation map source
        bounds : list
            The bounds of the object
        sid : int
            The source ID
        order : int
            The spectral order number
        """
        log.info(f"Dispersing source {sid}, order {order}")

        # Loop over all pixels in list for object "c"
        log.debug(f"source {sid} contains {len(xs)} pixels")
        all_res = []
        for i in range(len(xs)):
            # Here "i" indexes the pixel list for the object
            # being processed, as opposed to the ID number of the object itself

            width = 1.0
            height = 1.0
            xc = xs[i] + 0.5 * width
            yc = ys[i] + 0.5 * height

            fluxes, lams = map(np.array, zip(*[
                (fluxes_dict[lm][c][i], lm) for lm in sorted(fluxes_dict.keys())
                if fluxes_dict[lm][c][i] != 0
            ]))

            pars_i = (xc, yc, width, height, lams, fluxes, order,
                      wmin, wmax,
                      seg_wcs, grism_wcs, i, dims[::-1], 2,
                      extrapolate_sed, xoffset, yoffset)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in scalar divide")
                warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")
                all_res.append(dispersed_pixel(*pars_i))

        # Initialize blank image for this source
        this_object = np.zeros(dims, float)
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
            a = sparse.coo_matrix((f, (y - miny, x - minx)),
                                  shape=(maxy - miny + 1, maxx - minx + 1)).toarray()
            bounds.append([minx, maxx, miny, maxy])

            # Accumulate results into simulated images
            this_object[miny:maxy + 1, minx:maxx + 1] += a

            #if self.cache:
            #    self.cached_object[c]['x'].append(x)
            #    self.cached_object[c]['y'].append(y)
            #    self.cached_object[c]['f'].append(f)
            #    self.cached_object[c]['w'].append(w)
            #    self.cached_object[c]['minx'].append(minx)
            #    self.cached_object[c]['maxx'].append(maxx)
            #    self.cached_object[c]['miny'].append(miny)
            #    self.cached_object[c]['maxy'].append(maxy)

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
    def construct_slitmodel_for_chunk(chunk_data, bounds, sid, order):
        '''
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
        '''
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
        slit.data = chunk_data[thisobj_miny:thisobj_maxy + 1, thisobj_minx:thisobj_maxx + 1]

        return slit


    def disperse_all_from_cache(self, trans=None):
        if not self.cache:
            return

        self.simulated_image = np.zeros(self.dims, float)

        for i in range(len(self.IDs)):
            this_object = self.disperse_chunk_from_cache(i, trans=trans)

        return this_object

    def disperse_chunk_from_cache(self, c, trans=None):
        """Method that handles the dispersion. To be called after create_pixel_list()"""

        if not self.cache:
            return

        time1 = time.time()

        # Initialize blank image for this object
        this_object = np.zeros(self.dims, float)

        if trans is not None:
            log.debug("Applying a transmission function...")

        for i in range(len(self.cached_object[c]['x'])):
            x = self.cached_object[c]['x'][i]
            y = self.cached_object[c]['y'][i]
            f = self.cached_object[c]['f'][i] * 1.
            w = self.cached_object[c]['w'][i]

            if trans is not None:
                f *= trans(w)

            minx = self.cached_object[c]['minx'][i]
            maxx = self.cached_object[c]['maxx'][i]
            miny = self.cached_object[c]['miny'][i]
            maxy = self.cached_object[c]['maxy'][i]

            a = sparse.coo_matrix((f, (y - miny, x - minx)),
                                  shape=(maxy - miny + 1, maxx - minx + 1)).toarray()

            # Accumulate the results into the simulated images
            self.simulated_image[miny:maxy + 1, minx:maxx + 1] += a
            this_object[miny:maxy + 1, minx:maxx + 1] += a

        time2 = time.time()
        log.debug(f"Elapsed time {time2-time1} sec")

        return this_object
