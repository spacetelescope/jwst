"""
JWST pipeline step for image alignment.

:Authors: Mihai Cara

"""
from os import path
import math

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from tweakwcs.imalign import align_wcs
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.matchutils import XYXYMatch

from jwst.datamodels import ModelContainer

# LOCAL
from ..stpipe import Step
from ..assign_wcs.util import update_fits_wcsinfo, update_s_region_imaging, wcs_from_footprints
from .astrometric_utils import create_astrometric_catalog
from .tweakreg_catalog import make_tweakreg_catalog


_SQRT2 = math.sqrt(2.0)


def _oxford_or_str_join(str_list):
    nelem = len(str_list)
    if not nelem:
        return 'N/A'
    str_list = list(map(repr, str_list))
    if nelem == 1:
        return str_list
    elif nelem == 2:
        return str_list[0] + ' or ' + str_list[1]
    else:
        return ', '.join(map(repr, str_list[:-1])) + ', or ' + repr(str_list[-1])


SINGLE_GROUP_REFCAT = ['GAIADR3', 'GAIADR2', 'GAIADR1']
_SINGLE_GROUP_REFCAT_STR = _oxford_or_str_join(SINGLE_GROUP_REFCAT)

__all__ = ['TweakRegStep']


class TweakRegStep(Step):
    """
    TweakRegStep: Image alignment based on catalogs of sources detected in
    input images.
    """

    class_alias = "tweakreg"

    spec = f"""
        save_catalogs = boolean(default=False) # Write out catalogs?
        use_custom_catalogs = boolean(default=False) # Use custom user-provided catalogs?
        catalog_format = string(default='ecsv') # Catalog output file format
        catfile = string(default='') # Name of the file with a list of custom user-provided catalogs
        starfinder = option('dao', 'iraf', 'segmentation', default='iraf') # Star finder to use.

        # general starfinder options
        snr_threshold = float(default=10.0) # SNR threshold above the bkg for star finder
        bkg_boxsize = integer(default=400) # The background mesh box size in pixels.

        # kwargs for DAOStarFinder and IRAFStarFinder, only used if starfinder is 'dao' or 'iraf'
        kernel_fwhm = float(default=2.5) # Gaussian kernel FWHM in pixels
        minsep_fwhm = float(default=0.0) # Minimum separation between detected objects in FWHM
        sigma_radius = float(default=1.5) # Truncation radius of the Gaussian kernel in units of sigma
        sharplo = float(default=0.2) # The lower bound on sharpness for object detection.
        sharphi = float(default=1.0) # The upper bound on sharpness for object detection.
        roundlo = float(default=-1.0) # The lower bound on roundness for object detection.
        roundhi = float(default=1.0) # The upper bound on roundness for object detection.
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
        abs_refcat = string(default='')  # Catalog file name or one of: {_SINGLE_GROUP_REFCAT_STR}, or None, or ''
        save_abs_catalog = boolean(default=False)  # Write out used absolute astrometric reference catalog as a separate product

        # Absolute catalog align wcs options
        abs_minobj = integer(default=15) # Minimum number of objects acceptable for matching when performing absolute astrometry
        abs_fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift')
        abs_nclip = integer(min=0, default=3) # Number of clipping iterations in fit when performing absolute astrometry
        abs_sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units when performing absolute astrometry

        # absolute catalog xyxymatch options
        abs_searchrad = float(default=6.0) # The search radius in arcsec for a match when performing absolute astrometry
        # We encourage setting this parameter to True. Otherwise, xoffset and yoffset will be set to zero.
        abs_use2dhist = boolean(default=True) # Use 2D histogram to find initial offset when performing absolute astrometry?
        abs_separation = float(default=1) # Minimum object separation in arcsec when performing absolute astrometry
        abs_tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec when performing absolute astrometry

        # SIP approximation options, should match assign_wcs
        sip_approx = boolean(default=True)  # enables SIP approximation for imaging modes.
        sip_max_pix_error = float(default=0.01)  # max err for SIP fit, forward.
        sip_degree = integer(max=6, default=None)  # degree for forward SIP fit, None to use best fit.
        sip_max_inv_pix_error = float(default=0.01)  # max err for SIP fit, inverse.
        sip_inv_degree = integer(max=6, default=None)  # degree for inverse SIP fit, None to use best fit.
        sip_npoints = integer(default=12)  #  number of points for SIP
        
        # stpipe general options
        output_use_model = boolean(default=True)  # When saving use `DataModel.meta.filename`
    """

    reference_file_types = []

    def process(self, input):
        images = ModelContainer(input)

        if self.separation <= _SQRT2 * self.tolerance:
            self.log.error(
                "Parameter 'separation' must be larger than 'tolerance' by at "
                "least a factor of sqrt(2) to avoid source confusion."
            )
            for model in images:
                model.meta.cal_step.tweakreg = "SKIPPED"
            self.log.warning("Skipping 'TweakRegStep' step.")
            return input

        if self.abs_separation <= _SQRT2 * self.abs_tolerance:
            self.log.error(
                "Parameter 'abs_separation' must be larger than 'abs_tolerance' "
                "by at least a factor of sqrt(2) to avoid source confusion."
            )
            for model in images:
                model.meta.cal_step.tweakreg = "SKIPPED"
            self.log.warning("Skipping 'TweakRegStep' step.")
            return input

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
            elif hasattr(images.meta, "asn_table") and getattr(images, "asn_file_path", None) is not None:
                catdict = {}
                asn_dir = path.dirname(images.asn_file_path)
                for member in images.meta.asn_table.products[0].members:
                    if hasattr(member, "tweakreg_catalog"):
                        if member.tweakreg_catalog is None or not member.tweakreg_catalog.strip():
                            catdict[member.expname] = None
                        else:
                            catdict[member.expname] = path.join(asn_dir, member.tweakreg_catalog)

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
                self.log.warning("At least two exposures are required for image "
                                 "alignment.")
                self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                self.skip = True
                for model in images:
                    model.meta.cal_step.tweakreg = "SKIPPED"
                return input

        # === start processing images ===

        # pre-allocate collectors (same length and order as images)
        correctors = [None] * len(images)

        # Build the catalog and corrector for each input images
        for (model_index, image_model) in enumerate(images):
            # now that the model is open, check it's metadata for a custom catalog
            # only if it's not listed in the catdict
            if use_custom_catalogs and image_model.meta.filename not in catdict:
                if (image_model.meta.tweakreg_catalog is not None and image_model.meta.tweakreg_catalog.strip()):
                    catdict[image_model.meta.filename] = image_model.meta.tweakreg_catalog
            if use_custom_catalogs and catdict.get(image_model.meta.filename, None) is not None:
                # FIXME this modifies the input_model
                image_model.meta.tweakreg_catalog = catdict[image_model.meta.filename]
                # use user-supplied catalog:
                self.log.info("Using user-provided input catalog "
                              f"'{image_model.meta.tweakreg_catalog}'")
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
            catalog = _filter_catalog_by_bounding_box(
                catalog,
                image_model.meta.wcs.bounding_box)

            # setting 'name' is important for tweakwcs logging
            if catalog.meta.get('name') is None:
                catalog.meta['name'] = path.splitext(image_model.meta.filename)[0].strip('_- ')

            # log results of source finding (or user catalog)
            filename = image_model.meta.filename
            nsources = len(catalog)
            if nsources == 0:
                self.log.warning('No sources found in {}.'.format(filename))
            else:
                self.log.info('Detected {} sources in {}.'
                              .format(len(catalog), filename))

            # save catalog (if requested)
            if save_catalog:
                # FIXME this modifies the input_model
                image_model.meta.tweakreg_catalog = self._write_catalog(catalog, filename)

            # construct the corrector since the model is open (and already has a group_id)
            correctors[model_index] = _construct_wcs_corrector(image_model, catalog)

        self.log.info('')
        self.log.info("Number of image groups to be aligned: {:d}."
                      .format(n_groups))

        # keep track of if 'local' alignment failed, even if this
        # fails, absolute alignment might be run (if so configured)
        local_align_failed = False

        # if we have >1 group of images, align them to each other
        if n_groups > 1:

            # align images:
            xyxymatch = XYXYMatch(
                searchrad=self.searchrad,
                separation=self.separation,
                use2dhist=self.use2dhist,
                tolerance=self.tolerance,
                xoffset=self.xoffset,
                yoffset=self.yoffset
            )

            try:
                align_wcs(
                    correctors,
                    refcat=None,
                    enforce_user_order=self.enforce_user_order,
                    expand_refcat=self.expand_refcat,
                    minobj=self.minobj,
                    match=xyxymatch,
                    fitgeom=self.fitgeometry,
                    nclip=self.nclip,
                    sigma=(self.sigma, 'rmse')
                )

            except ValueError as e:
                msg = e.args[0]
                if (msg == "Too few input images (or groups of images) with "
                        "non-empty catalogs."):
                    # we need at least two exposures to perform image alignment
                    self.log.warning(msg)
                    self.log.warning("At least two exposures are required for "
                                     "image alignment.")
                    self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    if not align_to_abs_refcat:
                        self.skip = True
                        return images
                    local_align_failed = True
                else:
                    raise e

            except RuntimeError as e:
                msg = e.args[0]
                if msg.startswith("Number of output coordinates exceeded allocation"):
                    # we need at least two exposures to perform image alignment
                    self.log.error(msg)
                    self.log.error("Multiple sources within specified tolerance "
                                   "matched to a single reference source. Try to "
                                   "adjust 'tolerance' and/or 'separation' parameters.")
                    self.log.warning("Skipping 'TweakRegStep'...")
                    self.skip = True
                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    return images
                else:
                    raise e

            if not local_align_failed and not self._is_wcs_correction_small(correctors):
                if align_to_abs_refcat:
                    self.log.warning("Skipping relative alignment (stage 1)...")
                else:
                    self.log.warning("Skipping 'TweakRegStep'...")
                    self.skip = True
                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    return images

        if align_to_abs_refcat:
            # now, align things to the reference catalog
            # this can occur after alignment between groups (only if >1 group)

            # Get catalog of GAIA sources for the field
            #
            # NOTE:  If desired, the pipeline can write out the reference
            #        catalog as a separate product with a name based on
            #        whatever convention is determined by the JWST Cal Working
            #        Group.

            if self.save_abs_catalog:
                if self.output_dir is None:
                    output_name = 'fit_{}_ref.ecsv'.format(self.abs_refcat.lower())
                else:
                    output_name = path.join(self.output_dir, 'fit_{}_ref.ecsv'.format(self.abs_refcat.lower()))
            else:
                output_name = None

            # initial shift to be used with absolute astrometry
            self.abs_xoffset = 0
            self.abs_yoffset = 0

            self.abs_refcat = self.abs_refcat.strip()
            gaia_cat_name = self.abs_refcat.upper()

            if gaia_cat_name in SINGLE_GROUP_REFCAT:
                ref_model = images[0]

                epoch = Time(ref_model.meta.observation.date).decimalyear

                # combine all aligned wcs to compute a new footprint to
                # filter the absolute catalog sources
                combined_wcs = wcs_from_footprints(
                    None,
                    refmodel=ref_model,
                    wcslist=[corrector.wcs for corrector in correctors],
                )

                ref_cat = create_astrometric_catalog(
                    None,
                    gaia_cat_name,
                    existing_wcs=combined_wcs,
                    output=output_name,
                    epoch=epoch,
                )

            elif path.isfile(self.abs_refcat):
                ref_cat = Table.read(self.abs_refcat)

            else:
                raise ValueError("'abs_refcat' must be a path to an "
                                 "existing file name or one of the supported "
                                 f"reference catalogs: {_SINGLE_GROUP_REFCAT_STR}.")

            # Check that there are enough GAIA sources for a reliable/valid fit
            num_ref = len(ref_cat)
            if num_ref < self.abs_minobj:
                self.log.warning(
                    f"Not enough sources ({num_ref}) in the reference catalog "
                    "for the single-group alignment step to perform a fit. "
                    f"Skipping alignment to the {self.abs_refcat} reference "
                    "catalog!"
                )
            else:
                # align images:
                # Update to separation needed to prevent confusion of sources
                # from overlapping images where centering is not consistent or
                # for the possibility that errors still exist in relative overlap.
                xyxymatch_gaia = XYXYMatch(
                    searchrad=self.abs_searchrad,
                    separation=self.abs_separation,
                    use2dhist=self.abs_use2dhist,
                    tolerance=self.abs_tolerance,
                    xoffset=self.abs_xoffset,
                    yoffset=self.abs_yoffset
                )

                # Set group_id to same value so all get fit as one observation
                # The assigned value, 987654, has been hard-coded to make it
                # easy to recognize when alignment to GAIA was being performed
                # as opposed to the group_id values used for relative alignment
                # earlier in this step.
                for corrector in correctors:
                    corrector.meta['group_id'] = 987654
                    if ('fit_info' in corrector.meta and
                            'REFERENCE' in corrector.meta['fit_info']['status']):
                        del corrector.meta['fit_info']

                # Perform fit
                try:
                    align_wcs(
                        correctors,
                        refcat=ref_cat,
                        enforce_user_order=True,
                        expand_refcat=False,
                        minobj=self.abs_minobj,
                        match=xyxymatch_gaia,
                        fitgeom=self.abs_fitgeometry,
                        nclip=self.abs_nclip,
                        sigma=(self.abs_sigma, 'rmse')
                    )
                except ValueError as e:
                    msg = e.args[0]
                    if (msg == "Too few input images (or groups of images) with "
                            "non-empty catalogs."):
                        # we need at least two exposures to perform image alignment
                        self.log.warning(msg)
                        self.log.warning(
                            "At least one exposure is required to align images "
                            "to an absolute reference catalog. Alignment to an "
                            "absolute reference catalog will not be performed."
                        )
                        if local_align_failed or n_groups == 1:
                            self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
                            for model in images:
                                model.meta.cal_step.tweakreg = "SKIPPED"
                            self.skip = True
                            return images
                    else:
                        raise e

                except RuntimeError as e:
                    msg = e.args[0]
                    if msg.startswith("Number of output coordinates exceeded allocation"):
                        # we need at least two exposures to perform image alignment
                        self.log.error(msg)
                        self.log.error(
                            "Multiple sources within specified tolerance "
                            "matched to a single reference source. Try to "
                            "adjust 'tolerance' and/or 'separation' parameters."
                            "Alignment to an absolute reference catalog will "
                            "not be performed."
                        )
                        if local_align_failed or n_groups == 1:
                            self.log.warning("Skipping 'TweakRegStep'...")
                            self.skip = True
                            for model in images:
                                model.meta.cal_step.tweakreg = "SKIPPED"
                            return images
                    else:
                        raise e

        # one final pass through all the models to update them based
        # on the results of this step
        for (image_model, corrector) in zip(images, correctors):
            image_model.meta.cal_step.tweakreg = 'COMPLETE'

            # retrieve fit status and update wcs if fit is successful:
            if ('fit_info' in corrector.meta and
                    'SUCCESS' in corrector.meta['fit_info']['status']):

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
                    corrector.wcs.name = "FIT-LVL3-{}".format(self.abs_refcat)

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
                            crpix=None
                        )
                    except (ValueError, RuntimeError) as e:
                        self.log.warning(
                            "Failed to update 'meta.wcsinfo' with FITS SIP "
                            f'approximation. Reported error is:\n"{e.args[0]}"'
                        )

        return images

    def _write_catalog(self, catalog, filename):
        '''
        Determine output filename for catalog based on outfile for step
        and output dir, then write catalog to file.

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
        '''

        catalog_filename = str(filename).replace(
                    '.fits', '_cat.{}'.format(self.catalog_format)
                )
        if self.catalog_format == 'ecsv':
            fmt = 'ascii.ecsv'
        elif self.catalog_format == 'fits':
            # NOTE: The catalog must not contain any 'None' values.
            #       FITS will also not clobber existing files.
            fmt = 'fits'
        else:
            raise ValueError(
                '\'catalog_format\' must be "ecsv" or "fits".'
            )
        if self.output_dir is None:
            catalog.write(catalog_filename, format=fmt, overwrite=True)
        else:
            catalog.write(
                path.join(self.output_dir, catalog_filename),
                format=fmt,
                overwrite=True
            )
        self.log.info('Wrote source catalog: {}'
                      .format(catalog_filename))
        return catalog_filename

    def _find_sources(self, image_model):
        # source finding
        starfinder_kwargs = {
            'fwhm': self.kernel_fwhm,
            'sigma_radius': self.sigma_radius,
            'minsep_fwhm': self.minsep_fwhm,
            'sharplo': self.sharplo,
            'sharphi': self.sharphi,
            'roundlo': self.roundlo,
            'roundhi': self.roundhi,
            'peakmax': self.peakmax,
            'brightest': self.brightest,
            'npixels': self.npixels,
            'connectivity': int(self.connectivity),  # option returns a string, so cast to int
            'nlevels': self.nlevels,
            'contrast': self.contrast,
            'mode': self.multithresh_mode,
            'error': image_model.err,
            'localbkg_width': self.localbkg_width,
            'apermask_method': self.apermask_method,
            'kron_params': self.kron_params,
        }

        return make_tweakreg_catalog(
            image_model, self.snr_threshold,
            starfinder=self.starfinder,
            bkg_boxsize=self.bkg_boxsize,
            starfinder_kwargs=starfinder_kwargs,
        )

    def _is_wcs_correction_small(self, correctors):
        # check for a small wcs correction, it should be small
        if self.use2dhist:
            max_corr = 2 * (self.searchrad + self.tolerance) * u.arcsec
        else:
            max_corr = 2 * (max(abs(self.xoffset), abs(self.yoffset)) +
                            self.tolerance) * u.arcsec
        for corrector in correctors:
            aligned_skycoord = _wcs_to_skycoord(corrector.wcs)
            original_skycoord = corrector.meta['original_skycoord']
            separation = original_skycoord.separation(aligned_skycoord)
            if not (separation < max_corr).all():
                # Large corrections are typically a result of source
                # mis-matching or poorly-conditioned fit. Skip such models.
                self.log.warning(f"WCS has been tweaked by more than {10 * self.tolerance} arcsec")
                return False
        return True


def _parse_catfile(catfile):
    """
    Parse a text file containing at 2 whitespace-delimited columns
        column 1: str, datamodel filename
        column 2: str, catalog filename
    into a dictionary with datamodel filename keys and catalog filename
    values. The catalog filenames will become paths relative
    to the current working directory. So for a catalog filename
    "mycat.ecsv" if the catfile is in a subdirectory "my_data"
    the catalog filename will be "my_data/mycat.ecsv".

    Returns:
        - None of catfile is None (or an empty string)
        - empty dict if catfile is empty

    Raises:
        VaueError if catfile contains >2 columns
    """
    if catfile is None or not catfile.strip():
        return None

    catdict = {}

    with open(catfile) as f:
        catfile_dir = path.dirname(catfile)

        for line in f.readlines():
            sline = line.strip()
            if not sline or sline[0] == '#':
                continue

            data_model, *catalog = sline.split()
            catalog = list(map(str.strip, catalog))
            if len(catalog) == 1:
                catdict[data_model] = path.join(catfile_dir, catalog[0])
            elif len(catalog) == 0:
                # set this to None so it's custom catalog is skipped
                catdict[data_model] = None
            else:
                raise ValueError("'catfile' can contain at most two columns.")

    return catdict


def _rename_catalog_columns(catalog):
    for axis in ['x', 'y']:
        if axis not in catalog.colnames:
            long_axis = axis + 'centroid'
            if long_axis in catalog.colnames:
                catalog.rename_column(long_axis, axis)
            else:
                raise ValueError(
                    "'tweakreg' source catalogs must contain either "
                    "columns 'x' and 'y' or 'xcentroid' and "
                    "'ycentroid'."
                )
    return catalog


def _filter_catalog_by_bounding_box(catalog, bounding_box):
    if bounding_box is None:
        return catalog

    # filter out sources outside the WCS bounding box
    ((xmin, xmax), (ymin, ymax)) = bounding_box
    x = catalog['x']
    y = catalog['y']
    mask = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
    return catalog[mask]


def _wcs_to_skycoord(wcs):
    ra, dec = wcs.footprint(axis_type="spatial").T
    return SkyCoord(ra=ra, dec=dec, unit="deg")


def _construct_wcs_corrector(image_model, catalog):
    # pre-compute skycoord here so we can later use it
    # to check for a small wcs correction
    wcs = image_model.meta.wcs
    refang = image_model.meta.wcsinfo.instance
    return JWSTWCSCorrector(
        wcs=image_model.meta.wcs,
        wcsinfo={'roll_ref': refang['roll_ref'],
                 'v2_ref': refang['v2_ref'],
                 'v3_ref': refang['v3_ref']},
        # catalog and group_id are required meta
        meta={
            'catalog': catalog,
            'name': catalog.meta.get('name'),
            'group_id': image_model.meta.group_id,
            'original_skycoord': _wcs_to_skycoord(wcs),
        }
    )
