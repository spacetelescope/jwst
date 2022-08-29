"""
JWST pipeline step for image alignment.

:Authors: Mihai Cara

"""
from os import path

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from tweakwcs.imalign import align_wcs
from tweakwcs.correctors import JWSTWCSCorrector
from tweakwcs.matchutils import XYXYMatch

# LOCAL
from ..stpipe import Step
from .. import datamodels
from ..assign_wcs.util import update_fits_wcsinfo
from . import astrometric_utils as amutils
from .tweakreg_catalog import make_tweakreg_catalog


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


SINGLE_GROUP_REFCAT = ['GAIADR2', 'GAIADR1']
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
        catalog_format = string(default='ecsv') # Catalog output file format
        kernel_fwhm = float(default=2.5) # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=10.0) # SNR threshold above the bkg
        brightest = integer(default=200) # Keep top ``brightest`` objects
        peakmax = float(default=None) # Filter out objects with pixel values >= ``peakmax``
        enforce_user_order = boolean(default=False) # Align images in user specified order?
        expand_refcat = boolean(default=False) # Expand reference catalog with new sources?
        minobj = integer(default=15) # Minimum number of objects acceptable for matching
        searchrad = float(default=2.0) # The search radius in arcsec for a match
        use2dhist = boolean(default=True) # Use 2d histogram to find initial offset?
        separation = float(default=1.0) # Minimum object separation in arcsec
        tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec
        xoffset = float(default=0.0), # Initial guess for X offset in arcsec
        yoffset = float(default=0.0) # Initial guess for Y offset in arcsec
        fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift') # Fitting geometry
        nclip = integer(min=0, default=3) # Number of clipping iterations in fit
        sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units
        align_to_gaia = boolean(default=False)  # Align to GAIA catalog
        gaia_catalog = string(default='GAIADR2')  # Catalog file name or one of: {_SINGLE_GROUP_REFCAT_STR}
        min_gaia = integer(min=0, default=5) # Min number of GAIA sources needed
        save_gaia_catalog = boolean(default=False)  # Write out GAIA catalog as a separate product
        output_use_model = boolean(default=True)  # When saving use `DataModel.meta.filename`
        abs_minobj = integer(default=15) # Minimum number of objects acceptable for matching when performing absolute astrometry
        abs_searchrad = float(default=6.0) # The search radius in arcsec for a match when performing absolute astrometry
        # We encourage setting this parameter to True. Otherwise, xoffset and yoffset will be set to zero.
        abs_use2dhist = boolean(default=True) # Use 2D histogram to find initial offset when performing absolute astrometry?
        abs_separation = float(default=0.1) # Minimum object separation in arcsec when performing absolute astrometry
        abs_tolerance = float(default=0.7) # Matching tolerance for xyxymatch in arcsec when performing absolute astrometry
        # Fitting geometry when performing absolute astrometry
        abs_fitgeometry = option('shift', 'rshift', 'rscale', 'general', default='rshift')
        abs_nclip = integer(min=0, default=3) # Number of clipping iterations in fit when performing absolute astrometry
        abs_sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units when performing absolute astrometry
    """

    reference_file_types = []

    def process(self, input):

        try:
            images = datamodels.ModelContainer(input)
        except TypeError as e:
            e.args = ("Input to tweakreg must be a list of DataModels, an "
                      "association, or an already open ModelContainer "
                      "containing one or more DataModels.", ) + e.args[1:]
            raise e

        if self.align_to_gaia:
            # Set expand_refcat to True to eliminate possibility of duplicate
            # entries when aligning to GAIA
            self.expand_refcat = True

        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        # Build the catalogs for input images
        for image_model in images:
            # source finding
            catalog = make_tweakreg_catalog(
                image_model, self.kernel_fwhm, self.snr_threshold,
                brightest=self.brightest, peakmax=self.peakmax
            )

            # filter out sources outside the WCS bounding box
            bb = image_model.meta.wcs.bounding_box
            if bb is not None:
                ((xmin, xmax), (ymin, ymax)) = bb
                xname = 'xcentroid' if 'xcentroid' in catalog.colnames else 'x'
                yname = 'ycentroid' if 'ycentroid' in catalog.colnames else 'y'
                x = catalog[xname]
                y = catalog[yname]
                mask = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
                catalog = catalog[mask]

            filename = image_model.meta.filename
            nsources = len(catalog)
            if nsources == 0:
                self.log.warning('No sources found in {}.'.format(filename))
            else:
                self.log.info('Detected {} sources in {}.'
                              .format(len(catalog), filename))

            if self.save_catalogs:
                catalog_filename = filename.replace(
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
                catalog.write(catalog_filename, format=fmt, overwrite=True)
                self.log.info('Wrote source catalog: {}'
                              .format(catalog_filename))
                image_model.meta.tweakreg_catalog = catalog_filename

            # Temporarily attach catalog to the image model so that it follows
            # the grouping by exposure, to be removed after use below
            image_model.catalog = catalog

        # group images by their "group id":
        grp_img = list(images.models_grouped)

        self.log.info('')
        self.log.info("Number of image groups to be aligned: {:d}."
                      .format(len(grp_img)))
        self.log.info("Image groups:")

        if len(grp_img) == 1 and not self.align_to_gaia:
            self.log.info("* Images in GROUP 1:")
            for im in grp_img[0]:
                self.log.info("     {}".format(im.meta.filename))
            self.log.info('')

            # we need at least two exposures to perform image alignment
            self.log.warning("At least two exposures are required for image "
                             "alignment.")
            self.log.warning("Nothing to do. Skipping 'TweakRegStep'...")
            self.skip = True
            for model in images:
                model.meta.cal_step.tweakreg = "SKIPPED"
                # Remove the attached catalogs
                del model.catalog
            return input

        elif len(grp_img) > 1:

            # create a list of WCS-Catalog-Images Info and/or their Groups:
            imcats = []
            for g in grp_img:
                if len(g) == 0:
                    raise AssertionError("Logical error in the pipeline code.")
                else:
                    group_name = _common_name(g)
                    wcsimlist = list(map(self._imodel2wcsim, g))
                    # Remove the attached catalogs
                    for model in g:
                        del model.catalog
                    self.log.info("* Images in GROUP '{}':".format(group_name))
                    for im in wcsimlist:
                        im.meta['group_id'] = group_name
                        self.log.info("     {}".format(im.meta['name']))
                    imcats.extend(wcsimlist)

            self.log.info('')

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
                    imcats,
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
                    if not self.align_to_gaia:
                        self.skip = True
                        return images
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

            for imcat in imcats:
                model = imcat.meta['image_model']
                if model.meta.cal_step.tweakreg == "SKIPPED":
                    continue
                wcs = model.meta.wcs
                twcs = imcat.wcs
                if not self._is_wcs_correction_small(wcs, twcs):
                    # Large corrections are typically a result of source
                    # mis-matching or poorly-conditioned fit. Skip such models.
                    self.log.warning(f"WCS has been tweaked by more than {10 * self.tolerance} arcsec")

                    for model in images:
                        model.meta.cal_step.tweakreg = "SKIPPED"
                    if self.align_to_gaia:
                        self.log.warning("Skipping relative alignment (stage 1)...")
                    else:
                        self.log.warning("Skipping 'TweakRegStep'...")
                        self.skip = True
                        return images

        if self.align_to_gaia:
            # Get catalog of GAIA sources for the field
            #
            # NOTE:  If desired, the pipeline can write out the reference
            #        catalog as a separate product with a name based on
            #        whatever convention is determined by the JWST Cal Working
            #        Group.
            if self.save_gaia_catalog:
                output_name = 'fit_{}_ref.ecsv'.format(self.gaia_catalog.lower())
            else:
                output_name = None

            # initial shift to be used with absolute astrometry
            self.abs_xoffset = 0
            self.abs_yoffset = 0

            self.gaia_catalog = self.gaia_catalog.strip()
            gaia_cat_name = self.gaia_catalog.upper()

            if gaia_cat_name in SINGLE_GROUP_REFCAT:
                ref_cat = amutils.create_astrometric_catalog(
                    images,
                    gaia_cat_name,
                    output=output_name
                )

            elif path.isfile(self.gaia_catalog):
                ref_cat = Table.read(self.gaia_catalog)

            else:
                raise ValueError("'gaia_catalog' must be a path to an "
                                 "existing file name or one of the supported "
                                 f"reference catalogs: {_SINGLE_GROUP_REFCAT_STR}.")

            # Check that there are enough GAIA sources for a reliable/valid fit
            num_ref = len(ref_cat)
            if num_ref < self.min_gaia:
                # Raise Exception here to avoid rest of code in this try block
                self.log.warning(
                    f"Not enough sources ({num_ref}) in the reference catalog "
                    "for the single-group alignment step to perform a fit. "
                    f"Skipping alignment to the {self.gaia_catalog} reference "
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
                for imcat in imcats:
                    imcat.meta['group_id'] = 987654
                    if 'REFERENCE' in imcat.meta['fit_info']['status']:
                        del imcat.meta['fit_info']

                # Perform fit
                align_wcs(
                    imcats,
                    refcat=ref_cat,
                    enforce_user_order=True,
                    expand_refcat=False,
                    minobj=self.abs_minobj,
                    match=xyxymatch_gaia,
                    fitgeom=self.abs_fitgeometry,
                    nclip=self.abs_nclip,
                    sigma=(self.abs_sigma, 'rmse')
                )

        for imcat in imcats:
            image_model = imcat.meta['image_model']
            image_model.meta.cal_step.tweakreg = 'COMPLETE'

            # retrieve fit status and update wcs if fit is successful:
            if 'SUCCESS' in imcat.meta.get('fit_info')['status']:

                # Update/create the WCS .name attribute with information
                # on this astrometric fit as the only record that it was
                # successful:
                if self.align_to_gaia:
                    # NOTE: This .name attrib agreed upon by the JWST Cal
                    #       Working Group.
                    #       Current value is merely a place-holder based
                    #       on HST conventions. This value should also be
                    #       translated to the FITS WCSNAME keyword
                    #       IF that is what gets recorded in the archive
                    #       for end-user searches.
                    imcat.wcs.name = "FIT-LVL3-{}".format(self.gaia_catalog)

                image_model.meta.wcs = imcat.wcs

                # Also update FITS representation in input exposures for
                # subsequent reprocessing by the end-user.
                try:
                    update_fits_wcsinfo(
                        image_model,
                        max_pix_error=0.005
                    )
                except (ValueError, RuntimeError) as e:
                    self.log.warning(
                        "Failed to update 'meta.wcsinfo' with FITS SIP "
                        f'approximation. Reported error is:\n"{e.args[0]}"'
                    )

        return images

    def _is_wcs_correction_small(self, wcs, twcs):
        """Check that the newly tweaked wcs hasn't gone off the rails"""
        tolerance = 10.0 * self.tolerance * u.arcsec

        ra, dec = wcs.footprint(axis_type="spatial").T
        tra, tdec = twcs.footprint(axis_type="spatial").T
        skycoord = SkyCoord(ra=ra, dec=dec, unit="deg")
        tskycoord = SkyCoord(ra=tra, dec=tdec, unit="deg")

        separation = skycoord.separation(tskycoord)

        return (separation < tolerance).all()

    def _imodel2wcsim(self, image_model):
        # make sure that we have a catalog:
        if hasattr(image_model, 'catalog'):
            catalog = image_model.catalog
        else:
            catalog = image_model.meta.tweakreg_catalog

        model_name = path.splitext(image_model.meta.filename)[0].strip('_- ')

        if isinstance(catalog, Table):
            if not catalog.meta.get('name', None):
                catalog.meta['name'] = model_name

        else:
            try:
                cat_name = str(catalog)
                catalog = Table.read(catalog, format='ascii.ecsv')
                catalog.meta['name'] = cat_name
            except IOError:
                self.log.error("Cannot read catalog {}".format(catalog))

        if 'xcentroid' in catalog.colnames:
            catalog.rename_column('xcentroid', 'x')
            catalog.rename_column('ycentroid', 'y')

        # create WCSImageCatalog object:
        refang = image_model.meta.wcsinfo.instance
        im = JWSTWCSCorrector(
            wcs=image_model.meta.wcs,
            wcsinfo={'roll_ref': refang['roll_ref'],
                     'v2_ref': refang['v2_ref'],
                     'v3_ref': refang['v3_ref']},
            meta={'image_model': image_model, 'catalog': catalog,
                  'name': model_name}
        )

        return im


def _common_name(group):
    file_names = [path.splitext(im.meta.filename)[0].strip('_- ')
                  for im in group]
    fname_len = list(map(len, file_names))
    assert all(fname_len[0] == length for length in fname_len)
    cn = path.commonprefix(file_names)
    assert cn
    return cn
