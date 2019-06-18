"""
JWST pipeline step for image alignment.

:Authors: Mihai Cara

"""
from astropy.table import Table
from tweakwcs.imalign import align_wcs
from tweakwcs.tpwcs import JWSTgWCS
from tweakwcs.matchutils import TPMatch

# LOCAL
from ..stpipe import Step
from .. import datamodels

from .tweakreg_catalog import make_tweakreg_catalog


__all__ = ['TweakRegStep']


class TweakRegStep(Step):
    """
    TweakRegStep: Image alignment based on catalogs of sources detected in
    input images.
    """

    spec = """
        # Source finding parameters:
        save_catalogs = boolean(default=False) # Write out catalogs?
        catalog_format = string(default='ecsv') # Catalog output file format
        kernel_fwhm = float(default=2.5) # Gaussian kernel FWHM in pixels
        snr_threshold = float(default=10.0) # SNR threshold above the bkg
        brightest = integer(default=100) # Keep top ``brightest`` objects
        peakmax = float(default=None) # Filter out objects with pixel values >= ``peakmax``

        # Optimize alignment order:
        enforce_user_order = boolean(default=False) # Align images in user specified order?

        # Reference Catalog parameters:
        expand_refcat = boolean(default=False) # Expand reference catalog with new sources?

        # Object matching parameters:
        minobj = integer(default=15) # Minimum number of objects acceptable for matching
        searchrad = float(default=1.0) # The search radius in arcsec for a match
        use2dhist = boolean(default=True) # Use 2d histogram to find initial offset?
        separation = float(default=0.5) # Minimum object separation in arcsec
        tolerance = float(default=1.0) # Matching tolerance for xyxymatch in arcsec
        xoffset = float(default=0.0), # Initial guess for X offset in arcsec
        yoffset = float(default=0.0) # Initial guess for Y offset in arcsec

        # Catalog fitting parameters:
        fitgeometry = option('shift', 'rscale', 'general', default='general') # Fitting geometry
        nclip = integer(min=0, default=3) # Number of clipping iterations in fit
        sigma = float(min=0.0, default=3.0) # Clipping limit in sigma units
    """

    reference_file_types = []

    def process(self, input):

        try:
            images = datamodels.ModelContainer(input)
        except TypeError as te:
            raise te("Input to tweakreg must be a list of DataModels, an "
                "association, or an already open ModelContainer containing "
                "one or more DataModels.")

        # Build the catalogs for input images
        for image_model in images:
            catalog = make_tweakreg_catalog(
                image_model, self.kernel_fwhm, self.snr_threshold,
                brightest=self.brightest, peakmax=self.peakmax
            )

            # filter out sources outside the image array if WCS validity
            # region is provided:
            wcs_bounds = image_model.meta.wcs.pixel_bounds
            if wcs_bounds is not None:
                ((xmin, xmax), (ymin, ymax)) = wcs_bounds
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
                image_model.meta.tweakreg_catalog.filename = catalog_filename

            image_model.catalog = catalog

        # Now use the catalogs for tweakreg
        if len(images) == 0:
            raise ValueError("Input must contain at least one image model.")

        # group images by their "group id":
        grp_img = images.models_grouped

        if len(grp_img) == 1:
            # we need at least two exposures to perform image alignment
            self.log.info("At least two exposures are required for image "
                          "alignment.")
            self.log.info("Nothing to do. Skipping 'TweakRegStep'...")
            self.skip = True
            for model in images:
                model.meta.cal_step.tweakreg = "SKIPPED"
            return input

        # create a list of WCS-Catalog-Images Info and/or their Groups:
        imcats = []
        group_id = 1
        for g in grp_img:
            if len(g) == 0:
                raise AssertionError("Logical error in the pipeline code.")
            elif len(g) == 1:
                imcats.append(self._imodel2wcsim(g[0]))
            else:
                wcsimlist = list(map(self._imodel2wcsim, g))
                for im in wcsimlist:
                    im.meta['group_id'] = group_id
                imcats.extend(wcsimlist)
                group_id += 1

        # align images:
        tpmatch = TPMatch(
            searchrad=self.searchrad,
            separation=self.searchrad,
            use2dhist=self.use2dhist,
            tolerance=self.tolerance,
            xoffset=self.xoffset,
            yoffset=self.yoffset
        )

        align_wcs(
            imcats,
            refcat=None,
            enforce_user_order=self.enforce_user_order,
            expand_refcat=self.expand_refcat,
            minobj=self.minobj,
            match=tpmatch,
            fitgeom=self.fitgeometry,
            nclip=self.nclip,
            sigma=(self.sigma, 'rmse')
        )

        for imcat in imcats:
            fit_info = imcat.meta.get('fit_info')
            if fit_info['status'] == 'SUCCESS':
                imcat.meta['image_model'].meta.cal_step.tweakreg = 'COMPLETE'
                imcat.meta['image_model'].meta.wcs = imcat.wcs
            else:
                imcat.meta['image_model'].meta.cal_step.tweakreg = 'SKIPPED'

        return images

    def _imodel2wcsim(self, image_model):
        # make sure that we have a catalog:
        if hasattr(image_model, 'catalog'):
            catalog = image_model.catalog
        else:
            catalog = image_model.meta.tweakreg_catalog.filename

        if not isinstance(catalog, Table):
            try:
                catalog = Table.read(catalog, format='ascii.ecsv')
            except IOError:
                self.log.error("Cannot read catalog {}".format(catalog))

        if 'xcentroid' in catalog.colnames:
            catalog.rename_column('xcentroid', 'x')
            catalog.rename_column('ycentroid', 'y')

        # create WCSImageCatalog object:
        refang = image_model.meta.wcsinfo.instance
        im = JWSTgWCS(
            wcs=image_model.meta.wcs,
            wcsinfo={'roll_ref': refang['roll_ref'],
                     'v2_ref': refang['v2_ref'],
                     'v3_ref': refang['v3_ref']},
            meta={'image_model': image_model, 'catalog': catalog,
                  'name': image_model.meta.filename}
        )

        return im
