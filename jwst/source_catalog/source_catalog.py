import logging

import numpy as np

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip
from astropy.table import join, QTable
import astropy.units as u

from photutils import Background2D, MedianBackground
from photutils import detect_sources, deblend_sources, source_properties
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.utils import calc_total_error
from photutils.utils._wcs_helpers import _pixel_scale_angle_at_skycoord

from ..datamodels import ImageModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class SourceCatalog:
    def __init__(self, model, bkg_boxsize=100, kernel_fwhm=2.0,
                 snr_threshold=3.0, npixels=5.0, deblend=False,
                 aperture_ee=(30, 50, 70)):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')

        self.model = model
        self.bkg_boxsize = np.int(bkg_boxsize)  # must be integer
        self.kernel_fwhm = kernel_fwhm
        self.snr_threshold = snr_threshold
        self.npixels = npixels
        self.deblend = deblend
        self.aperture_ee = aperture_ee

        self._flux_colnames = self._make_aper_colnames('flux')
        self._abmag_colnames = self._make_aper_colnames('abmag')
        self._vegamag_colnames = self._make_aper_colnames('vegamag')

        self.coverage_mask = None
        self._bkg = None
        self.background = None
        self.background_rms = None
        self.kernel = None
        self.threshold = None
        self.segm = None
        self.total_error = None
        self.segment_catalog = None
        self.null_column = None

        self.bkg_aper_in = None
        self.bkg_aper_out = None
        self.aper_radii = None
        self.abmag_offset = None
        self.aperture_catalog = None
        self._ci_colname = None

    def _colnames_generator(self, name):
        for ee in self.aperture_ee:
            base = f'aper{ee}_{name}'
            yield base
            yield f'{base}_err'

    def _make_aper_colnames(self, name):
        colnames = list(self._colnames_generator(name))
        colnames.extend([f'aper_total_{name}', f'aper_total_{name}_err'])
        return colnames

    def convert_to_jy(self):
        """ convert from MJy/sr to Jy"""

        if self.model.meta.bunit_data != 'MJy/sr':
            raise ValueError('data is expected to be in units of MJy/sr')
        self.model.data *= (1.e6 *
                            self.model.meta.photometry.pixelarea_steradians)
        self.model.meta.bunit_data = 'Jy'

        return

    def make_coverage_mask(self):
        # NOTE: currently the wht output from the resample step is
        # always an exposure-time map
        self.coverage_mask = (self.model.wht == 0)

        return self.coverage_mask

    def estimate_background(self):
        """
        Estimate the 2D background and background RMS noise in an image.

        Parameters
        ----------
        box_size : int or array_like (int)
            The box size along each axis.  If ``box_size`` is a scalar then
            a square box of size ``box_size`` will be used.  If ``box_size``
            has two elements, they should be in ``(ny, nx)`` order.

        coverage_mask : bool ndarray

        Returns
        -------
        background : `photutils.background.Background2D`
            A Background2D object containing the 2D background and
            background RMS noise estimates.
        """

        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = MedianBackground()
        filter_size = (3, 3)
        bkg = Background2D(self.model.data, self.bkg_boxsize,
                           filter_size=filter_size, mask=self.coverage_mask,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        # apply the coverage mask
        bkg.background *= ~self.coverage_mask
        bkg.background_rms *= ~self.coverage_mask

        self._bkg = bkg
        self.background = bkg.background
        self.background_rms = bkg.background_rms

        return self.background, self.background_rms

    def subtract_background(self):
        self.model.data -= self.background
        return

    def make_kernel(self):
        """
        Make a 2D Gaussian smoothing kernel that is used to filter the image
        before thresholding.

        Filtering the image will smooth the noise and maximize detectability
        of objects with a shape similar to the kernel.

        Parameters
        ----------
        kernel_fwhm : float
            The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.

        Returns
        -------
        kernel : `astropy.convolution.Kernel2D`
            The output smoothing kernel, normalized such that it sums to 1.
        """

        sigma = self.kernel_fwhm * gaussian_fwhm_to_sigma
        #kernel = Gaussian2DKernel(sigma)
        kernel = Gaussian2DKernel(sigma, x_size=5, y_size=5)
        kernel.normalize(mode='integral')

        self.kernel = kernel

        return self.kernel

    def calc_threshold(self):
        self.threshold = self.snr_threshold * self.background_rms
        return self.threshold

    def detect_sources(self):
        """
        Detect sources in an image, including deblending.

        Parameters
        ----------
        model : `ImageModel`
            The input `ImageModel` of a single drizzled image.  The
            input image is assumed to be background subtracted.

        threshold : float
            The data value or pixel-wise data values to be used for the
            detection threshold. A 2D threshold must have the same shape as
            ``model.data``.

        npixels : int
            The number of connected pixels, each greater than the threshold
            that an object must have to be detected.  ``npixels`` must be a
            positive integer.

        kernel : `astropy.convolution.Kernel2D`
            The filtering kernel.

        mask : bool ndarray

        deblend : bool, optional
            Whether to deblend overlapping sources.  Source deblending
            requires scikit-image.

        Returns
        -------
        segment_image : `~photutils.segmentation.SegmentationImage` or `None`
            A 2D segmentation image, with the same shape as the input data
            where sources are marked by different positive integer values.
            A value of zero is reserved for the background.  If no sources
            are found then `None` is returned.
        """

        connectivity = 8
        segm = detect_sources(self.model.data, self.threshold,
                              npixels=self.npixels, filter_kernel=self.kernel,
                              mask=self.coverage_mask,
                              connectivity=connectivity)

        # segm=None for photutils >= 0.7
        # segm.nlabels=0 for photutils < 0.7
        if segm is None or segm.nlabels == 0:
            return None

        # source deblending requires scikit-image
        if self.deblend:
            nlevels = 32
            contrast = 0.001
            mode = 'exponential'
            segm = deblend_sources(self.model.data, segm,
                                   npixels=self.npixels,
                                   filter_kernel=self.kernel, nlevels=nlevels,
                                   contrast=contrast, mode=mode,
                                   connectivity=connectivity, relabel=True)
        self.segm = segm

        return segm

    def _calc_total_error(self):
        """
        Calculate total error (background plus source Poisson error).
        """

        # The resample step still does not allow for background-only
        # inverse-variance maps.  Here we estimate the background-only rms
        # directly from the data, but note this is dependent on the
        # ``box_size`` parameter in ``estimate_background``.
        bkg_error = self.background_rms  # Jy

        # The source Poisson noise contribution to the total error requires
        # converting the image into electrons.  The input images are in
        # units of MJy/sr.  They can be converted to units of counts (DN)
        # per second via the PHOTMJSR keyword, and then to counts (DN) using
        # the exposure-time weight map.  However, the final conversion from
        # DN to electrons requires knowing the gain, which is not available
        # in the level-3 products.  Here we assume the gain is 1.

        # the effective gain is the multiplicative factor to convert the
        # input image from Jy to e-
        gain = 1.0  # e-/DN; this is a rough approximation
        dn_s_factor = (1.e6 *
                       self.model.meta.photometry.pixelarea_steradians *
                       self.model.meta.photometry.conversion_megajanskys)
        effective_gain = gain / dn_s_factor * self.model.wht  # Jy -> e-

        self.total_error = calc_total_error(self.model.data, bkg_error,
                                            effective_gain)

        return self.total_error

    def calc_total_error(self):
        self.total_error = np.zeros_like(self.model.data)
        return self.total_error

    def apply_units(self):
        unit = u.Jy
        self.model.data <<= unit
        self.total_error <<= unit
        self.background <<= unit
        self.background_rms <<= unit
        return

    @staticmethod
    def flux_to_abmag(flux, flux_err):
        """
        Parameters
        ----------
        flux : float, ndarray
            In units of Jy
        """

        abmag = -2.5 * np.log10(flux.value) + 8.9
        # assuming SNR >> 1 (otherwise abmag_err is asymmetric)
        abmag_err = 2.5 * np.log10(np.e) * (flux_err.value / flux.value)

        return abmag, abmag_err

    def make_segment_catalog(self):
        """
        Create a catalog of source photometry and morphologies based on
        segmentation.
        """

        # wcs = self.model.meta.wcs  # gWCS
        wcs = self.model.get_fits_wcs()  # FITS WCS

        source_props = source_properties(self.model.data.astype(float),
                                         self.segm, error=self.total_error,
                                         background=self.background,
                                         filter_kernel=self.kernel, wcs=wcs)

        columns = ['id', 'xcentroid', 'ycentroid', 'sky_centroid',
                   'source_sum', 'source_sum_err', 'area',
                   'semimajor_axis_sigma', 'semiminor_axis_sigma',
                   'ellipticity', 'orientation', 'sky_bbox_ll', 'sky_bbox_ul',
                   'sky_bbox_lr', 'sky_bbox_ur']
        catalog = source_props.to_table(columns=columns)

        flux_col = 'isophotal_flux'
        flux_err_col = 'isophotal_flux_err'
        catalog.rename_column('source_sum', flux_col)
        catalog.rename_column('source_sum_err', flux_err_col)
        catalog.rename_column('area', 'isophotal_area')
        catalog.rename_column('semimajor_axis_sigma', 'semimajor_sigma')
        catalog.rename_column('semiminor_axis_sigma', 'semiminor_sigma')

        # define the orientation position angle (E of N)
        # NOTE: crpix1 and crpix2 are 1-based values
        skycoord = wcs.pixel_to_world(self.model.meta.wcsinfo.crpix1 - 1,
                                      self.model.meta.wcsinfo.crpix2 - 1)
        _, angle = _pixel_scale_angle_at_skycoord(skycoord, wcs)
        sky_orientation = (180.0 * u.deg) - angle + catalog['orientation']
        catalog.add_column(sky_orientation, name='sky_orientation', index=11)

        flux = catalog[flux_col]
        flux_err = catalog[flux_err_col]
        abmag, abmag_err = self.flux_to_abmag(flux, flux_err)
        vegamag = abmag - self.abmag_offset
        vegamag_err = abmag_err
        catalog.add_column(abmag, name='isophotal_abmag', index=6)
        catalog.add_column(abmag_err, name='isophotal_abmag_err', index=7)
        catalog.add_column(vegamag, name='isophotal_vegamag', index=8)
        catalog.add_column(vegamag_err, name='isophotal_vegamag_err', index=9)

        self.segment_catalog = catalog
        self.xypos = np.transpose((catalog['xcentroid'],
                                   catalog['ycentroid']))

        return self.segment_catalog

    def make_null_column(self):
        values = np.empty(len(self.segment_catalog))
        values.fill(np.nan)
        self.null_column = values

    def calc_aper_local_background(self, bkg_aper):
        """
        Estimate the local background as the sigma-clipped median value
        in the input aperture.
        """

        bkg_aper_masks = bkg_aper.to_mask(method='center')
        sigclip = SigmaClip(sigma=3)

        bkg_median = []
        bkg_std = []
        nvalues = []
        for mask in bkg_aper_masks:
            bkg_data = mask.multiply(self.model.data.value)
            bkg_data_1d = bkg_data[mask.data > 0]
            values = sigclip(bkg_data_1d, masked=False)
            nvalues.append(values.size)
            bkg_median.append(np.median(values))
            bkg_std.append(np.std(values))

        nvalues = np.array(nvalues)
        bkg_median = np.array(bkg_median)
        # standard error of the median
        bkg_median_err = np.sqrt(np.pi / (2. * nvalues)) * np.array(bkg_std)

        bkg_median <<= self.model.data.unit
        bkg_median_err <<= self.model.data.unit

        return bkg_median, bkg_median_err

    def append_aper_total(self, catalog):
        catalog['aper_total_flux'] = self.null_column
        catalog['aper_total_flux_err'] = self.null_column
        catalog['aper_total_abmag'] = self.null_column
        catalog['aper_total_abmag_err'] = self.null_column
        catalog['aper_total_vegamag'] = self.null_column
        catalog['aper_total_vegamag_err'] = self.null_column

        return catalog

    def calc_concentration_index(self, catalog):
        eefractions = [self.aperture_ee[i] for i in (0, 2)]
        self._ci_colname = f'CI_{eefractions[0]}_{eefractions[1]}'
        col1 = f'aper{eefractions[0]}_abmag'
        col2 = f'aper{eefractions[1]}_abmag'
        return catalog[col1] - catalog[col2]

    def calc_is_star(self, catalog):
        return self.null_column

    def get_aperture_radii(self):
        self.bkg_aper_in = 10
        self.bkg_aper_out = 15
        self.aper_radii = (3, 5, 7)
        return

    def get_abmag_offset(self):
        filepath = 'abmag_to_vega.ecsv'
        tbl = QTable.read(filepath)
        self.abmag_offset = 0.
        return self.abmag_offset

    def make_aperture_catalog(self):
        self.get_aperture_radii()
        bkg_aperture = CircularAnnulus(self.xypos, self.bkg_aper_in,
                                       self.bkg_aper_out)
        bkg_median, bkg_median_err = self.calc_aper_local_background(
            bkg_aperture)

        apertures = [CircularAperture(self.xypos, radius) for radius in
                     self.aper_radii]
        aper_phot = aperture_photometry(self.model.data, apertures,
                                        error=self.total_error)

        for i, aperture in enumerate(apertures):
            flux_col = f'aperture_sum_{i}'
            flux_err_col = f'aperture_sum_err_{i}'

            # subtract the local background measured in the annulus
            aper_phot[flux_col] -= (bkg_median * aperture.area)

            flux = aper_phot[flux_col]
            flux_err = aper_phot[flux_err_col]
            abmag, abmag_err = self.flux_to_abmag(flux, flux_err)
            vegamag = abmag - self.abmag_offset
            vegamag_err = abmag_err

            idx0 = 2 * i
            idx1 = (2 * i) + 1
            new_flux_col = self._flux_colnames[idx0]
            new_flux_err_col = self._flux_colnames[idx1]
            abmag_col = self._abmag_colnames[idx0]
            abmag_err_col = self._abmag_colnames[idx1]
            vegamag_col = self._vegamag_colnames[idx0]
            vegamag_err_col = self._vegamag_colnames[idx1]

            aper_phot.rename_column(flux_col, new_flux_col)
            aper_phot.rename_column(flux_err_col, new_flux_err_col)

            aper_phot[abmag_col] = abmag
            aper_phot[abmag_err_col] = abmag_err
            aper_phot[vegamag_col] = vegamag
            aper_phot[vegamag_err_col] = vegamag_err

        aper_phot['aper_bkg_flux'] = bkg_median
        aper_phot['aper_bkg_fluxerr'] = bkg_median_err

        aper_phot = self.append_aper_total(aper_phot)

        aper_phot[self._ci_colname] = self.calc_concentration_index(aper_phot)
        aper_phot['is_star'] = self.calc_is_star(aper_phot)

        self.aperture_catalog = aper_phot

        return self.aperture_catalog

    def calc_sharpness(self, catalog):
        return self.null_column

    def calc_roundness(self, catalog):
        return self.null_column

    def calc_isolation_metric(self, catalog):
        return self.null_column

    def make_source_catalog(self):
        catalog = join(self.segment_catalog, self.aperture_catalog)

        catalog['sharpness'] = self.calc_sharpness(catalog)
        catalog['roundness'] = self.calc_roundness(catalog)
        catalog['isolation_metric'] = self.calc_isolation_metric(catalog)

        for colname in catalog.colnames:
            if colname.endswith('_err'):
                catalog[colname] = self.null_column

        colnames = self.segment_catalog.colnames[0:4]
        colnames.extend(self._flux_colnames)
        colnames.extend(['aper_bkg_flux', 'aper_bkg_fluxerr'])
        colnames.extend(self._abmag_colnames)
        colnames.extend(self._vegamag_colnames)
        colnames.extend([self._ci_colname, 'is_star', 'sharpness',
                         'roundness', 'isolation_metric'])
        colnames.extend(self.segment_catalog.colnames[4:])

        self.catalog = catalog[colnames]

        return self.catalog

    def run(self):
        self.convert_to_jy()
        self.make_coverage_mask()
        self.estimate_background()
        self.subtract_background()
        self.make_kernel()
        self.calc_threshold()

        segm = self.detect_sources()
        if segm is None:
            log.info('No sources were found. Source catalog will not be '
                     'created.')
            return None
        log.info(f'Detected {segm.nlabels} sources')

        self.calc_total_error()
        self.apply_units()
        self.get_abmag_offset()

        self.make_segment_catalog()
        self.make_nan_column()
        self.make_aperture_catalog()
        self.make_source_catalog()

        return self.catalog
