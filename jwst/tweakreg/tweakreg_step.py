"""
JWST pipeline step for image alignment.

:Authors: Mihai Cara
"""

from pathlib import Path

from astropy.table import Table
from astropy.time import Time
from tweakwcs.correctors import JWSTWCSCorrector

import stcal.tweakreg.tweakreg as twk

from jwst.stpipe import record_step_status
from jwst.assign_wcs.util import update_fits_wcsinfo, update_s_region_imaging
from jwst.datamodels import ModelLibrary

# LOCAL
from jwst.stpipe import Step
from .tweakreg_catalog import make_tweakreg_catalog


def _oxford_or_str_join(str_list):
    nelem = len(str_list)
    if not nelem:
        return "N/A"
    str_list = list(map(repr, str_list))
    if nelem == 1:
        return str_list
    elif nelem == 2:
        return str_list[0] + " or " + str_list[1]
    else:
        return ", ".join(map(repr, str_list[:-1])) + ", or " + repr(str_list[-1])


SINGLE_GROUP_REFCAT = ["GAIADR3", "GAIADR2", "GAIADR1"]
_SINGLE_GROUP_REFCAT_STR = _oxford_or_str_join(SINGLE_GROUP_REFCAT)

__all__ = ["TweakRegStep"]


class TweakRegStep(Step):
    """Image alignment based on catalogs of sources detected in input images."""

    class_alias = "tweakreg"

    spec = f"""
        save_catalogs = boolean(default=False) # Write out catalogs?
        use_custom_catalogs = boolean(default=False) # Use custom user-provided catalogs?
        catalog_format = string(default='ecsv') # Catalog output file format
        catfile = string(default='') # Name of the file with a list of custom user-provided catalogs
        starfinder = option('dao', 'iraf', 'segmentation', default='iraf') # Star finder to use.

        # general starfinder options
        snr_threshold = float(default=10.0) # SNR threshold above the bkg for star finder
        kernel_fwhm = float(default=2.5) # Gaussian kernel FWHM in pixels
        bkg_boxsize = integer(default=400) # The background mesh box size in pixels.

        # kwargs for DAOStarFinder and IRAFStarFinder, only used if starfinder is 'dao' or 'iraf'
        minsep_fwhm = float(default=0.0) # Minimum separation between detected objects in FWHM
        sigma_radius = float(default=1.5) # Truncation radius of the Gaussian kernel, units of sigma
        sharplo = float(default=0.5) # The lower bound on sharpness for object detection.
        sharphi = float(default=2.0) # The upper bound on sharpness for object detection.
        roundlo = float(default=0.0) # The lower bound on roundness for object detection.
        roundhi = float(default=0.2) # The upper bound on roundness for object detection.
        brightest = integer(default=200) # Keep top ``brightest`` objects
        peakmax = float(default=None) # Filter out objects with pixel values >= ``peakmax``

        # kwargs for SourceCatalog and SourceFinder, only used if starfinder is 'segmentation'
        npixels = integer(default=10) # Minimum number of connected pixels
        connectivity = option(4, 8, default=8) # The connectivity defining the neighborhood of a pixel
        nlevels = integer(default=32) # Number of multi-thresholding levels for deblending
        contrast = float(default=0.001) # Fraction of total source flux an object must have to be deblended
        multithresh_mode = option('exponential', 'linear', 'sinh', default='exponential') # Multi-thresholding mode
        localbkg_width = integer(default=0) # Width of rectangular annulus used to compute local background around each source
        apermask_method = option('correct', 'mask', 'none', default='correct') # How to handle neighboring sources
        kron_params = float_list(min=2, max=3, default=None) # Parameters defining Kron aperture

        # align wcs options
        enforce_user_order = boolean(default=False) # Align images in user specified order?
        expand_refcat = boolean(default=False) # Expand reference catalog with new sources?
        minobj = integer(default=15) # Minimum number of objects acceptable for matching
        fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift') # Fitting geometry
        nclip = integer(min=0, default=3) # Number of clipping iterations in fit
        sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units

        # xyxymatch options
        searchrad = float(default=2.0) # The search radius in arcsec for a match
        use2dhist = boolean(default=True) # Use 2d histogram to find initial offset?
        separation = float(default=1.0) # Minimum object separation for xyxymatch in arcsec
        tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec
        xoffset = float(default=0.0), # Initial guess for X offset in arcsec
        yoffset = float(default=0.0) # Initial guess for Y offset in arcsec

        # Absolute catalog options
        abs_refcat = string(default='') # Catalog file name or one of: {_SINGLE_GROUP_REFCAT_STR}, or None, or ''
        save_abs_catalog = boolean(default=False) # Write out used absolute astrometric reference catalog as a separate product

        # Absolute catalog align wcs options
        abs_minobj = integer(default=15) # Minimum number of objects acceptable for matching when performing absolute astrometry
        abs_fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift')
        abs_nclip = integer(min=0, default=3) # Number of clipping iterations in fit when performing absolute astrometry
        abs_sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units when performing absolute astrometry

        # absolute catalog xyxymatch options
        abs_searchrad = float(default=6.0) # The search radius in arcsec for a match when performing absolute astrometry
        abs_use2dhist = boolean(default=True) # Use 2D histogram to find initial offset when performing absolute astrometry? We encourage setting this parameter to True. Otherwise, xoffset and yoffset will be set to zero.
        abs_separation = float(default=1) # Minimum object separation in arcsec when performing absolute astrometry
        abs_tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec when performing absolute astrometry

        # SIP approximation options, should match assign_wcs
        sip_approx = boolean(default=True) # enables SIP approximation for imaging modes.
        sip_max_pix_error = float(default=0.01) # max err for SIP fit, forward.
        sip_degree = integer(max=6, default=None) # degree for forward SIP fit, None to use best fit.
        sip_max_inv_pix_error = float(default=0.01) # max err for SIP fit, inverse.
        sip_inv_degree = integer(max=6, default=None) # degree for inverse SIP fit, None to use best fit.
        sip_npoints = integer(default=12) #  number of points for SIP

        # stpipe general options
        output_use_model = boolean(default=True) # When saving use `DataModel.meta.filename`
        in_memory = boolean(default=True) # If False, preserve memory using temporary files at expense of runtime
    """  # noqa: E501

    reference_file_types: list = []

    def process(self, input_data):
        """
        Perform image alignment based on catalogs of sources detected in input images.

        Parameters
        ----------
        input_data : `ModelLibrary`, or asn-type input to be read into a `ModelLibrary`
            A collection of data models.

        Returns
        -------
        output : `ModelLibrary`
            The aligned input data models.
        """
        if isinstance(input_data, ModelLibrary):
            images = input_data
        else:
            images = ModelLibrary(input_data, on_disk=not self.in_memory)

        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        # determine number of groups (used below)
        n_groups = len(images.group_names)

        use_custom_catalogs = self.use_custom_catalogs

        if self.use_custom_catalogs:
            # first check catfile
            if self.catfile.strip():
                catdict = _parse_catfile(self.catfile)
                # if user requested the use of custom catalogs and provided a
                # valid 'catfile' file name that has no custom catalogs,
                # turn off the use of custom catalogs:
                if not catdict:
                    self.log.warning(
                        "'use_custom_catalogs' is set to True but 'catfile' "
                        "contains no user catalogs. Turning on built-in catalog "
                        "creation."
                    )
                    use_custom_catalogs = False
            # else, load from association
            elif images.asn_dir is not None:
                catdict = {}
                for member in images.asn["products"][0]["members"]:
                    if "tweakreg_catalog" in member:
                        tweakreg_catalog = member["tweakreg_catalog"]
                        if tweakreg_catalog is None or not tweakreg_catalog.strip():
                            catdict[member["expname"]] = None
                        else:
                            # convert back to string to allow schema to validate
                            cat = str(Path(images.asn_dir) / tweakreg_catalog)
                            catdict[member["expname"]] = cat

        if self.abs_refcat is not None and self.abs_refcat.strip():
            align_to_abs_refcat = True
            # Set expand_refcat to True to eliminate possibility of duplicate
            # entries when aligning to absolute astrometric reference catalog
            self.expand_refcat = True
        else:
            align_to_abs_refcat = False

            # since we're not aligning to a reference catalog, check if we
            # are saving catalogs, if not, and we have 1 group, skip
            if not self.save_catalogs and n_groups == 1:
                # we need at least two exposures to perform image alignment
                self.log.warning("At least two exposures are required for image alignment.")
                self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                record_step_status(images, "tweakreg", success=False)
                return images

        # === start processing images ===

        # pre-allocate collectors (same length and order as images)
        correctors = [None] * len(images)

        # Build the catalog and corrector for each input images
        with images:
            for model_index, image_model in enumerate(images):
                # now that the model is open, check its metadata for a custom catalog
                # only if it's not listed in the catdict
                if use_custom_catalogs and image_model.meta.filename not in catdict:
                    if (
                        image_model.meta.tweakreg_catalog is not None
                        and image_model.meta.tweakreg_catalog.strip()
                    ):
                        catdict[image_model.meta.filename] = image_model.meta.tweakreg_catalog
                if use_custom_catalogs and catdict.get(image_model.meta.filename, None) is not None:
                    image_model.meta.tweakreg_catalog = catdict[image_model.meta.filename]
                    # use user-supplied catalog:
                    self.log.info(
                        f"Using user-provided input catalog '{image_model.meta.tweakreg_catalog}'"
                    )
                    catalog = Table.read(
                        image_model.meta.tweakreg_catalog,
                    )
                    save_catalog = False
                else:
                    # source finding
                    catalog = self._find_sources(image_model)

                    # only save if catalog was computed from _find_sources and
                    # the user requested save_catalogs
                    save_catalog = self.save_catalogs

                # if needed rename xcentroid to x, ycentroid to y
                catalog = _rename_catalog_columns(catalog)

                # filter all sources outside the wcs bounding box
                catalog = twk.filter_catalog_by_bounding_box(
                    catalog, image_model.meta.wcs.bounding_box
                )

                # setting 'name' is important for tweakwcs logging
                if catalog.meta.get("name") is None:
                    catalog.meta["name"] = Path(image_model.meta.filename).stem.strip("_- ")

                # log results of source finding (or user catalog)
                filename = image_model.meta.filename
                nsources = len(catalog)
                if nsources == 0:
                    self.log.warning(f"No sources found in {filename}.")
                else:
                    self.log.info(f"Detected {len(catalog)} sources in {filename}.")

                # save catalog (if requested)
                if save_catalog:
                    # FIXME this modifies the input_model
                    image_model.meta.tweakreg_catalog = self._write_catalog(catalog, filename)

                # construct the corrector since the model is open (and already has a group_id)
                correctors[model_index] = twk.construct_wcs_corrector(
                    image_model.meta.wcs,
                    image_model.meta.wcsinfo.instance,
                    catalog,
                    image_model.meta.group_id,
                )
                images.shelve(image_model, model_index)

        self.log.info("")
        self.log.info(f"Number of image groups to be aligned: {n_groups}.")

        # wrapper to stcal tweakreg routines
        # step skip conditions should throw TweakregError from stcal
        if n_groups > 1:
            try:
                # relative alignment of images to each other (if more than one group)
                correctors = twk.relative_align(
                    correctors,
                    enforce_user_order=self.enforce_user_order,
                    expand_refcat=self.expand_refcat,
                    minobj=self.minobj,
                    fitgeometry=self.fitgeometry,
                    nclip=self.nclip,
                    sigma=self.sigma,
                    searchrad=self.searchrad,
                    use2dhist=self.use2dhist,
                    separation=self.separation,
                    tolerance=self.tolerance,
                    xoffset=self.xoffset,
                    yoffset=self.yoffset,
                )
            except twk.TweakregError as e:
                self.log.warning(str(e))
                local_align_failed = True
            else:
                local_align_failed = False
        else:
            local_align_failed = True

        # absolute alignment to the reference catalog
        # can (and does) occur after alignment between groups
        if align_to_abs_refcat:
            self.log.info(f"Aligning to absolute reference catalog: {self.abs_refcat}")
            with images:
                ref_image = images.borrow(0)
                try:
                    correctors = twk.absolute_align(
                        correctors,
                        self.abs_refcat,
                        ref_wcs=ref_image.meta.wcs,
                        ref_wcsinfo=ref_image.meta.wcsinfo.instance,
                        epoch=Time(ref_image.meta.observation.date).decimalyear,
                        abs_minobj=self.abs_minobj,
                        abs_fitgeometry=self.abs_fitgeometry,
                        abs_nclip=self.abs_nclip,
                        abs_sigma=self.abs_sigma,
                        abs_searchrad=self.abs_searchrad,
                        abs_use2dhist=self.abs_use2dhist,
                        abs_separation=self.abs_separation,
                        abs_tolerance=self.abs_tolerance,
                        save_abs_catalog=self.save_abs_catalog,
                        abs_catalog_output_dir=self.output_dir,
                    )
                    images.shelve(ref_image, 0, modify=False)
                except twk.TweakregError as e:
                    self.log.warning(str(e))
                    images.shelve(ref_image, 0, modify=False)
                    record_step_status(images, "tweakreg", success=False)
                    return images
                finally:
                    del ref_image

        if local_align_failed and not align_to_abs_refcat:
            record_step_status(images, "tweakreg", success=False)
            return images

        # one final pass through all the models to update them based
        # on the results of this step
        self._apply_tweakreg_solution(images, correctors, align_to_abs_refcat=align_to_abs_refcat)
        return images

    def _apply_tweakreg_solution(
        self,
        images: ModelLibrary,
        correctors: list[JWSTWCSCorrector],
        align_to_abs_refcat: bool = False,
    ) -> ModelLibrary:
        """
        Apply the WCS corrections to the input images.

        Parameters
        ----------
        images : ModelLibrary
            A collection of data models.
        correctors : list[JWSTWCSCorrector]
            A list of WCS correctors.
        align_to_abs_refcat : bool
            Flag indicating whether the images were aligned to an absolute reference catalog.

        Returns
        -------
        images : ModelLibrary
            The aligned input data models
        """
        with images:
            for image_model, corrector in zip(images, correctors, strict=True):
                # retrieve fit status and update wcs if fit is successful:
                if (
                    "fit_info" in corrector.meta
                    and "SUCCESS" in corrector.meta["fit_info"]["status"]
                ):
                    # Update/create the WCS .name attribute with information
                    # on this astrometric fit as the only record that it was
                    # successful:
                    if align_to_abs_refcat:
                        # NOTE: This .name attrib agreed upon by the JWST Cal
                        #       Working Group.
                        #       Current value is merely a place-holder based
                        #       on HST conventions. This value should also be
                        #       translated to the FITS WCSNAME keyword
                        #       IF that is what gets recorded in the archive
                        #       for end-user searches.
                        corrector.wcs.name = f"FIT-LVL3-{self.abs_refcat}"

                    image_model.meta.wcs = corrector.wcs
                    update_s_region_imaging(image_model)

                    # Also update FITS representation in input exposures for
                    # subsequent reprocessing by the end-user.
                    if self.sip_approx:
                        try:
                            update_fits_wcsinfo(
                                image_model,
                                max_pix_error=self.sip_max_pix_error,
                                degree=self.sip_degree,
                                max_inv_pix_error=self.sip_max_inv_pix_error,
                                inv_degree=self.sip_inv_degree,
                                npoints=self.sip_npoints,
                                crpix=None,
                            )
                        except (ValueError, RuntimeError) as e:
                            self.log.warning(
                                "Failed to update 'meta.wcsinfo' with FITS SIP "
                                "approximation. Reported error is:"
                            )
                            self.log.warning(f'"{e.args[0]}"')
                record_step_status(image_model, "tweakreg", success=True)
                images.shelve(image_model)
        return images

    def _write_catalog(self, catalog, filename):
        """
        Determine output filename and write catalog to file.

        Parameters
        ----------
        catalog : astropy.table.Table
            Table containing the source catalog.
        filename : str
            Output filename for step

        Returns
        -------
        catalog_filename : str
            Filename where the catalog was saved
        """
        catalog_filename = str(filename).replace(".fits", f"_cat.{self.catalog_format}")
        if self.catalog_format == "ecsv":
            fmt = "ascii.ecsv"
        elif self.catalog_format == "fits":
            # NOTE: The catalog must not contain any 'None' values.
            #       FITS will also not clobber existing files.
            fmt = "fits"
        else:
            raise ValueError('\'catalog_format\' must be "ecsv" or "fits".')
        if self.output_dir is None:
            catalog.write(catalog_filename, format=fmt, overwrite=True)
        else:
            catalog.write(Path(self.output_dir) / catalog_filename, format=fmt, overwrite=True)
        self.log.info(f"Wrote source catalog: {catalog_filename}")
        return catalog_filename

    def _find_sources(self, image_model):
        # source finding
        starfinder_kwargs = {
            "sigma_radius": self.sigma_radius,
            "minsep_fwhm": self.minsep_fwhm,
            "sharplo": self.sharplo,
            "sharphi": self.sharphi,
            "roundlo": self.roundlo,
            "roundhi": self.roundhi,
            "peakmax": self.peakmax,
            "brightest": self.brightest,
            "npixels": self.npixels,
            "connectivity": int(self.connectivity),  # option returns a string, so cast to int
            "nlevels": self.nlevels,
            "contrast": self.contrast,
            "mode": self.multithresh_mode,
            "error": image_model.err,
            "localbkg_width": self.localbkg_width,
            "apermask_method": self.apermask_method,
            "kron_params": self.kron_params,
        }

        return make_tweakreg_catalog(
            image_model,
            self.snr_threshold,
            self.kernel_fwhm,
            starfinder_name=self.starfinder,
            bkg_boxsize=self.bkg_boxsize,
            starfinder_kwargs=starfinder_kwargs,
        )


def _parse_catfile(catfile):
    """
    Parse catalog text file into a dictionary.

    The catalog text file must contain two whitespace-delimited columns:
        column 1: str, datamodel filename
        column 2: str, catalog filename
    The catalog filenames will become paths relative
    to the current working directory. So for a catalog filename
    "mycat.ecsv" if the catfile is in a subdirectory "my_data"
    the catalog filename will be "my_data/mycat.ecsv".

    Parameters
    ----------
    catfile : str
        Path to a text file containing the list of catalogs.

    Returns
    -------
    catdict : dict or None
        Dictionary with datamodel filename keys and catalog filename values.
        None if catfile is None (or an empty string).
        empty dict if catfile is empty.

    Raises
    ------
    ValueError if catfile contains >2 columns
    """
    if catfile is None or not catfile.strip():
        return None

    catdict = {}

    catfile = Path(catfile)
    with catfile.open() as f:
        catfile_dir = catfile.parent

        for line in f.readlines():
            sline = line.strip()
            if not sline or sline[0] == "#":
                continue

            data_model, *catalog = sline.split()
            catalog = list(map(str.strip, catalog))
            if len(catalog) == 1:
                # convert back to string to allow schema to validate
                catdict[data_model] = str(catfile_dir / catalog[0])
            elif len(catalog) == 0:
                # set this to None so it's custom catalog is skipped
                catdict[data_model] = None
            else:
                raise ValueError("'catfile' can contain at most two columns.")

    return catdict


def _rename_catalog_columns(catalog):
    """
    Rename columns 'xcentroid' and 'ycentroid' to 'x' and 'y', respectively.

    Parameters
    ----------
    catalog : astropy.table.Table
        Table containing the source catalog.

    Returns
    -------
    catalog : astropy.table.Table
        Table containing the source catalog with renamed columns.
    """
    for axis in ["x", "y"]:
        if axis not in catalog.colnames:
            long_axis = axis + "centroid"
            if long_axis in catalog.colnames:
                catalog.rename_column(long_axis, axis)
            else:
                raise ValueError(
                    "'tweakreg' source catalogs must contain either "
                    "columns 'x' and 'y' or 'xcentroid' and "
                    "'ycentroid'."
                )
    return catalog
