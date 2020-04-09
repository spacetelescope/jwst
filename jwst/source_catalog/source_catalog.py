import logging

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip
from astropy.table import join
import astropy.units as u
import numpy as np
from scipy.spatial import cKDTree

from photutils import Background2D, MedianBackground
from photutils import detect_sources, deblend_sources, source_properties
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.utils import calc_total_error
from photutils.utils._wcs_helpers import _pixel_scale_angle_at_skycoord

from .. import datamodels
from ..datamodels import ImageModel, ABVegaOffsetModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class SourceCatalog:
    def __init__(self, model, bkg_boxsize=100, kernel_fwhm=2.0,
                 snr_threshold=3.0, npixels=5.0, deblend=False,
                 aperture_ee=(30, 50, 70), apcorr_filename=None,
                 abvega_offset_filename=None):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')

        self.model = model
        self.bkg_boxsize = np.int(bkg_boxsize)  # must be integer
        self.kernel_fwhm = kernel_fwhm
        self.snr_threshold = snr_threshold
        self.npixels = npixels
        self.deblend = deblend
        self.apcorr_filename = apcorr_filename
        self.abvega_offset_filename = abvega_offset_filename

        aperture_ee = np.array(aperture_ee).astype(int)
        if len(aperture_ee) != 3:
            raise ValueError('aperture_ee must contain only 3 values')
        if np.any(np.logical_or(aperture_ee <= 0, aperture_ee >= 100)):
            raise ValueError('aperture_ee values must be between 0 and 100')
        if np.any(aperture_ee[1:] < aperture_ee[:-1]):
            raise ValueError('aperture_ee values must be in increasing order')
        self.aperture_ee = aperture_ee

        self.instrument = self.model.meta.instrument.name
        self.detector = self.model.meta.instrument.detector
        self.filtername = self.model.meta.instrument.filter
        self.pupil = model.meta.instrument.pupil
        self.subarray = self.model.meta.subarray.name
        log.info(f'Instrument: {self.instrument}')
        if self.detector is not None:
            log.info(f'Detector: {self.detector}')
        if self.filtername is not None:
            log.info(f'Filter: {self.filtername}')
        if self.pupil is not None:
            log.info(f'Pupil: {self.pupil}')
        if self.subarray is not None:
            log.info(f'Subarray: {self.subarray}')

        self._ci_colname = f'CI_{self.aperture_ee[0]}_{self.aperture_ee[2]}'
        self._flux_colnames, self._flux_coldesc = self._make_aperture_colnames('flux')
        self._abmag_colnames, self._abmag_coldesc = self._make_aperture_colnames('abmag')
        self._vegamag_colnames, self._vegamag_coldesc = self._make_aperture_colnames('vegamag')
        self._bkg_colnames = ['aper_bkg_flux', 'aper_bkg_flux_err']
        self._bkg_coldesc = ['The local background value calculated as the sigma-clipped median value in the background annulus aperture', 'The standard error of the sigma-clipped median background value']

        self.coverage_mask = None
        self._bkg = None
        self.background = None
        self.background_rms = None
        self.kernel = None
        self.threshold = None
        self.segm = None
        self.total_error = None
        self.segment_catalog = None
        self.xypos = None
        self.null_column = None

        self.is_star = None
        self.aperture_radii = None
        self.aperture_corrs = None
        self.bkg_aperture_inner = None
        self.bkg_aperture_outer = None
        self.abvega_offset = None
        self.aperture_catalog = None
        self.extras_catalog = None
        self.catalog = None

    @staticmethod
    def _add_column(catalog, column, value, description=None):
        catalog[column] = value
        if description is not None:
            catalog[column].info.description = description
        return

    def _make_aperture_colnames(self, name):
        if name == 'flux':
            ftype = 'Flux'
            ftype2 = 'flux'
        elif name == 'abmag':
            ftype = ftype2 = 'AB magnitude'
        elif name == 'vegamag':
            ftype = ftype2 = 'Vega magnitude'

        colnames = []
        descriptions = []
        for aper_ee in self.aperture_ee:
            basename = f'aper{aper_ee}_{name}'
            colnames.append(basename)
            descriptions.append(f'{ftype} within the {aper_ee}% encircled energy circular aperture')
            colnames.append(f'{basename}_err')
            descriptions.append(f'{ftype} error within the {aper_ee}% encircled energy circular aperture')

        colnames.extend([f'aper_total_{name}', f'aper_total_{name}_err'])
        descriptions.append(f'Total aperture-corrected {ftype2} based on the {self.aperture_ee[-1]}% encircled energy circular aperture - calculated only for stars')
        descriptions.append(f'Total aperture-corrected {ftype2} error based on the {self.aperture_ee[-1]}% encircled energy circular aperture - calculated only for stars')

        return colnames, descriptions

    def _find_aperture_params(self, selector):
        apcorr_model = datamodels.open(self.apcorr_filename)
        apcorr = apcorr_model.apcorr_table
        if selector is None:  # FGS
            ee_table = apcorr
        else:
            mask_idx = [apcorr[key] == value
                        for key, value in selector.items()]
            ee_table = apcorr[np.logical_and.reduce(mask_idx)]

        radii = []
        apcorrs = []
        skyins = []
        skyouts = []
        for ee in self.aperture_ee:
            row = ee_table[np.round(ee_table['eefraction'] * 100) == ee]
            if len(row) == 0:
                raise RuntimeError('Aperture encircled energy value of {0} '
                                   'appears to be invalid. No matching row '
                                   'was found in the apcorr reference file '
                                   '{1}'.format(ee, self.apcorr_filename))
            if len(row) > 1:
                raise RuntimeError('More than one matching row was found in '
                                   'the apcorr reference file {0}'
                                   .format(self.apcorr_filename))

            radii.append(row['radius'][0])
            apcorrs.append(row['apcorr'][0])
            skyins.append(row['skyin'][0])
            skyouts.append(row['skyout'][0])

        # TODO
        # These aperture and background-aperture radii are in pixel
        # units and are based on the detector native pixel scales.  If
        # the user changes the output pixel scale in the resample step,
        # then these radii need to be scaled.  Changing the output pixel
        # scale is not yet a pipeline option (as of Apr 2020).
        self.aperture_radii = np.array(radii)
        self.aperture_corrs = np.array(apcorrs)

        skyins = np.unique(skyins)
        skyouts = np.unique(skyouts)
        if len(skyins) != 1 or len(skyouts) != 1:
            raise RuntimeError('Expected to find only one value for skyin '
                               'and skyout in the apcorr reference file for '
                               'a given selector.')
        self.bkg_aperture_inner = skyins[0]
        self.bkg_aperture_outer = skyouts[0]

        return

    def set_aperture_params(self):
        if self.instrument == 'NIRCAM' or self.instrument == 'NIRISS':
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername, 'subarray': self.subarray}
        elif self.instrument == 'FGS':
            selector = None
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        self._find_aperture_params(selector)

        return

    def set_abvega_offset(self):
        if self.instrument == 'NIRCAM' or self.instrument == 'NIRISS':
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername}
        elif self.instrument == 'FGS':
            selector = {'detector': self.detector}
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        if self.abvega_offset_filename is None:
            self.abvega_offset = 0.0
            log.info('ABVegaOffsetModel reference file was not input -- '
                     'catalog Vega magnitudes are not correct.')
            return

        abvega_offset_model = ABVegaOffsetModel(self.abvega_offset_filename)
        offsets_table = abvega_offset_model.abvega_offset

        try:
            mask_idx = [offsets_table[key] == value
                        for key, value in selector.items()]
        except KeyError as badkey:
            raise KeyError('{0} not found in ABVegaOffsetModel reference '
                           'file {1}'.format(badkey,
                                             self.abvega_offset_filename))

        row = offsets_table[np.logical_and.reduce(mask_idx)]

        if len(row) == 0:
            raise RuntimeError('Did not find matching row in '
                               'ABVegaOffsetModel reference file {0}'
                               .format(self.abvega_offset_filename))
        if len(row) > 1:
            raise RuntimeError('Found more than one matching row in '
                               'ABVegaOffsetModel reference file {0}'
                               .format(self.abvega_offset_filename))

        self.abvega_offset = row['abvega_offset'][0]
        log.info('AB to Vega magnitude offset {:.5f}'
                 .format(self.abvega_offset))

        return

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
        bkg.background *= np.logical_not(self.coverage_mask)
        bkg.background_rms *= np.logical_not(self.coverage_mask)

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
        kernel = Gaussian2DKernel(sigma)
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
        #self.total_error = self._calc_total_error()
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

        # handle negative fluxes
        idx = flux.value < 0
        abmag[idx] = np.nan
        abmag_err[idx] = np.nan

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

        # save the source positions for aperture photometry
        self.xypos = np.transpose((catalog['xcentroid'],
                                   catalog['ycentroid']))

        flux_col = 'isophotal_flux'
        flux_err_col = 'isophotal_flux_err'
        catalog.rename_column('source_sum', flux_col)
        catalog.rename_column('source_sum_err', flux_err_col)
        catalog.rename_column('area', 'isophotal_area')
        catalog.rename_column('semimajor_axis_sigma', 'semimajor_sigma')
        catalog.rename_column('semiminor_axis_sigma', 'semiminor_sigma')

        # Define the orientation position angle (E of N)
        # NOTE: crpix1 and crpix2 are 1-based values
        skycoord = wcs.pixel_to_world(self.model.meta.wcsinfo.crpix1 - 1,
                                      self.model.meta.wcsinfo.crpix2 - 1)
        _, angle = _pixel_scale_angle_at_skycoord(skycoord, wcs)
        sky_orientation = (180.0 * u.deg) - angle + catalog['orientation']
        catalog.add_column(sky_orientation, name='sky_orientation', index=11)

        # Add ABmag and Vegamag columns
        flux = catalog[flux_col]
        flux_err = catalog[flux_err_col]
        abmag, abmag_err = self.flux_to_abmag(flux, flux_err)
        vegamag = abmag - self.abvega_offset
        vegamag_err = abmag_err
        catalog.add_column(abmag, name='isophotal_abmag', index=6)
        catalog.add_column(abmag_err, name='isophotal_abmag_err', index=7)
        catalog.add_column(vegamag, name='isophotal_vegamag', index=8)
        catalog.add_column(vegamag_err, name='isophotal_vegamag_err', index=9)

        # add column descriptions
        desc = {}
        desc['id'] = 'Unique source identification number'
        desc['xcentroid'] = 'X pixel value of the source centroid'
        desc['ycentroid'] = 'Y pixel value of the source centroid'
        desc['sky_centroid'] = 'Sky coordinate of the source centroid'
        desc['isophotal_flux'] = 'Isophotal flux'
        desc['isophotal_flux_err'] = 'Isophotal flux error'
        desc['isophotal_abmag'] = 'Isophotal AB magnitude'
        desc['isophotal_abmag_err'] = 'Isophotal AB magnitude error'
        desc['isophotal_vegamag'] = 'Isophotal Vega magnitude'
        desc['isophotal_vegamag_err'] = 'Isophotal Vega magnitude error'
        desc['isophotal_area'] = 'Isophotal area'
        desc['semimajor_sigma'] = ('1-sigma standard deviation along the '
                                   'semimajor axis of the 2D Gaussian '
                                   'function that has the same second-order '
                                   'central moments as the source')
        desc['semiminor_sigma'] = ('1-sigma standard deviation along the '
                                   'semiminor axis of the 2D Gaussian '
                                   'function that has the same second-order '
                                   'central moments as the source')
        desc['ellipticity'] = ('1 minus the ratio of the 1-sigma lengths of '
                               'the semimajor and semiminor axes')
        desc['orientation'] = ('The angle (degrees) between the positive X '
                               'axis and the major axis (increases '
                               'counter-clockwise)')
        desc['sky_orientation'] = ('The position angle (degrees) from North '
                                   'of the major axis')
        desc['sky_bbox_ll'] = ('Sky coordinate of the lower-left vertex of '
                               'the minimal bounding box of the source')
        desc['sky_bbox_ul'] = ('Sky coordinate of the upper-left vertex of '
                               'the minimal bounding box of the source')
        desc['sky_bbox_lr'] = ('Sky coordinate of the lower-right vertex of '
                               'the minimal bounding box of the source')
        desc['sky_bbox_ur'] = ('Sky coordinate of the upper-right vertex of '
                               'the minimal bounding box of the source')

        for key, val in desc.items():
            catalog[key].info.description = val


        #catalog['id'].info.description = 'Unique source identification number'
        #catalog['xcentroid'].info.description = 'X pixel value of the source centroid'
        #catalog['ycentroid'].info.description = 'Y pixel value of the source centroid'
        #catalog['sky_centroid'].info.description = 'Sky coordinate of the source centroid'
        #catalog['isophotal_flux'].info.description = 'Isophotal flux'
        #catalog['isophotal_flux_err'].info.description = 'Isophotal flux error'
        #catalog['isophotal_abmag'].info.description = 'Isophotal AB magnitude'
        #catalog['isophotal_abmag_err'].info.description = 'Isophotal AB magnitude error'
        #catalog['isophotal_vegamag'].info.description = 'Isophotal Vega magnitude'
        #catalog['isophotal_vegamag_err'].info.description = 'Isophotal Vega magnitude error'
        #catalog['isophotal_area'].info.description = 'Isophotal area'
        #catalog['semimajor_sigma'].info.description = '1-sigma standard deviation along the semimajor axis of the 2D Gaussian function that has the same second-order central moments as the source'
        #catalog['semiminor_sigma'].info.description = '1-sigma standard deviation along the semiminor axis of the 2D Gaussian function that has the same second-order central moments as the source'
        #catalog['ellipticity'].info.description = '1 minus the ratio of the 1-sigma lengths of the semimajor and semiminor axes'
        #catalog['orientation'].info.description = 'The angle (degrees) between the positive X axis and the major axis (increases counter-clockwise)'
        #catalog['sky_orientation'].info.description = 'The position angle (degrees) from North of the major axis'
        #catalog['sky_bbox_ll'].info.description = 'Sky coordinate of the lower-left vertex of the minimal bounding box of the source'
        #catalog['sky_bbox_ul'].info.description = 'Sky coordinate of the upper-left vertex of the minimal bounding box of the source'
        #catalog['sky_bbox_lr'].info.description = 'Sky coordinate of the lower-right vertex of the minimal bounding box of the source'
        #catalog['sky_bbox_ur'].info.description = 'Sky coordinate of the upper-right vertex of the minimal bounding box of the source'

        self.segment_catalog = catalog

        return self.segment_catalog

    def make_null_column(self):
        values = np.empty(len(self.segment_catalog))
        values.fill(np.nan)
        self.null_column = values

        return self.null_column

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

    def calc_concentration_index(self):
        col1 = self._abmag_colnames[0]  # ee_fraction[0]
        col2 = self._abmag_colnames[4]  # ee_fraction[2]

        return self.aperture_catalog[col1] - self.aperture_catalog[col2]

    def calc_is_star(self):
        """
        Calculate a flag based on whether the source is a star.

        2020.03.02: The JWST Photometry Working Group has not determined
        the criteria for this flag.
        """

        # TODO: need algorithm for this flag
        #is_star = np.random.randint(2, size=len(self.aperture_catalog))
        #return is_star.astype(bool)

        return self.null_column

    def append_aper_total(self, catalog):
        """
        Calculate the aperture-corrected total flux for sources that are
        stars.

        The aperture fluxes in the largest aperture (i.e., the largest
        encircled energy) are corrected to total flux.
        """

        idx = 2  # apcorr for the largest EE (largest radius)
        apcorr = self.aperture_corrs[idx]
        total_flux = (apcorr *
                      self.aperture_catalog[self._flux_colnames[idx*2]])
        total_flux_err = (apcorr *
                          self.aperture_catalog[self._flux_colnames[idx*2 + 1]])

        if np.isnan(self.is_star[0]):
            # TODO: remove this when ``is_star`` values are available
            # this array is all False
            is_star = np.zeros(len(self.aperture_catalog), dtype=bool)
        else:
            is_star = self.is_star

        not_star = np.logical_not(is_star)
        total_flux[not_star] = np.nan
        total_flux_err[not_star] = np.nan
        catalog[self._flux_colnames[-2]] = total_flux
        catalog[self._flux_colnames[-1]] = total_flux_err

        abmag, abmag_err = self.flux_to_abmag(total_flux, total_flux_err)
        vegamag = abmag - self.abvega_offset
        vegamag_err = abmag_err
        #catalog[self._abmag_colnames[-2]] = abmag
        #catalog[self._abmag_colnames[-1]] = abmag_err
        #catalog[self._vegamag_colnames[-2]] = vegamag
        #catalog[self._vegamag_colnames[-1]] = vegamag_err

        for idx in (-2, -1):
            self._add_column(catalog, self._flux_colnames[idx], total_flux,
                             self._flux_coldesc[idx])
            self._add_column(catalog, self._abmag_colnames[idx], total_flux,
                             self._abmag_coldesc[idx])
            self._add_column(catalog, self._vegamag_colnames[idx], total_flux,
                             self._vegamag_coldesc[idx])


        #catalog[self._flux_colnames[-2]].info.description = self._flux_coldesc[-2]
        #catalog[self._flux_colnames[-1]].info.description = self._flux_coldesc[-1]
        #catalog[self._abmag_colnames[-2]].info.description = self._abmag_coldesc[-2]
        #catalog[self._abmag_colnames[-1]].info.description = self._abmag_coldesc[-1]
        #catalog[self._vegamag_colnames[-2]].info.description = self._vegamag_coldesc[-2]
        #catalog[self._vegamag_colnames[-1]].info.description = self._vegamag_coldesc[-2]

        return catalog

    def make_aperture_catalog(self):
        bkg_aperture = CircularAnnulus(self.xypos, self.bkg_aperture_inner,
                                       self.bkg_aperture_outer)
        bkg_median, bkg_median_err = self.calc_aper_local_background(
            bkg_aperture)

        apertures = [CircularAperture(self.xypos, radius) for radius in
                     self.aperture_radii]
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
            vegamag = abmag - self.abvega_offset
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
            aper_phot[new_flux_col].info.description = self._flux_coldesc[idx0]
            aper_phot[new_flux_err_col].info.description = self._flux_coldesc[idx1]

            #print(abmag_col, abmag, self._abmag_coldesc[idx0])

            self._add_column(aper_phot, abmag_col, abmag,
                             self._abmag_coldesc[idx0])
            self._add_column(aper_phot, abmag_err_col, abmag_err,
                             self._abmag_coldesc[idx1])
            self._add_column(aper_phot, vegamag_col, vegamag,
                             self._vegamag_coldesc[idx0])
            self._add_column(aper_phot, vegamag_err_col, vegamag_err,
                             self._vegamag_coldesc[idx1])

        aper_phot.add_column(bkg_median, name=self._bkg_colnames[0], index=3)
        aper_phot.add_column(bkg_median_err, name=self._bkg_colnames[1],
                             index=4)
        aper_phot[self._bkg_colnames[0]].info.description = self._bkg_coldesc[0]
        aper_phot[self._bkg_colnames[1]].info.description = self._bkg_coldesc[1]

        self.aperture_catalog = aper_phot

        return self.aperture_catalog

    def calc_sharpness(self, catalog):
        return self.null_column

    def calc_roundness(self, catalog):
        return self.null_column

    def calc_nearest_neighbors(self, catalog):
        """
        Calculate the distance in pixels to the nearest neighbor and
        its AB magnitude.

        If the nearest-neighbor is a star, then its aperture-corrected
        total AB magnitude is used.  Otherwise, its isophotal AB
        magnitude is used.
        """

        xypos = np.transpose((self.segment_catalog['xcentroid'],
                              self.segment_catalog['ycentroid']))
        tree = cKDTree(xypos)
        qdist, qidx = tree.query(xypos, k=2)
        dist = np.transpose(qdist)[1]
        idx = np.transpose(qidx)[1]

        nn_abmag = catalog['aper_total_abmag'][idx]
        nn_isomag = self.segment_catalog['isophotal_abmag'][idx]

        return dist, np.where(np.isnan(nn_abmag), nn_isomag, nn_abmag)

    def make_extras_catalog(self):
        catalog = self.segment_catalog[['id']]  # new table

        self.is_star = self.calc_is_star()
        catalog = self.append_aper_total(catalog)

        ci_desc = ('Concentration index calculated as '
                   f'{self._abmag_colnames[0]} - {self._abmag_colnames[4]}')
        self._add_column(catalog, self._ci_colname,
                         self.calc_concentration_index(), ci_desc)

        is_star_desc = 'Flag indicating whether the source is a star'
        self._add_column(catalog, 'is_star', self.is_star, is_star_desc)

        sharpness_desc = 'The DAOFind source sharpness statistic'
        self._add_column(catalog, 'sharpness', self.calc_sharpness(catalog),
                         sharpness_desc)

        roundness_desc = 'The DAOFind source roundness statistic'
        self._add_column(catalog, 'roundness', self.calc_roundness(catalog),
                         roundness_desc)

        nn_dist, nn_abmag = self.calc_nearest_neighbors(catalog)
        nn_dist_desc = 'The distance in pixels to the nearest neighbor'
        self._add_column(catalog, 'nn_dist', nn_dist * u.pixel, nn_dist_desc)

        nn_abmag_desc = 'The AB magnitude of the nearest neighbor'
        self._add_column(catalog, 'nn_abmag', nn_abmag, nn_abmag_desc)

        self.extras_catalog = catalog

        return self.extras_catalog

    def make_source_catalog(self):
        catalog = join(self.segment_catalog, self.aperture_catalog)
        catalog = join(catalog, self.extras_catalog)

        # TODO: Until errors are produced for the level 3 drizzle
        # products, the JPWG has decided that the errors should not be
        # populated.
        for colname in catalog.colnames:
            if colname.endswith('_err'):
                catalog[colname][:] = self.null_column

        # Define the column name order for the final source catalog
        colnames = self.segment_catalog.colnames[0:4]
        colnames.extend(self._bkg_colnames)
        colnames.extend(self._flux_colnames)
        colnames.extend(self._abmag_colnames)
        colnames.extend(self._vegamag_colnames)
        colnames.extend(self.extras_catalog.colnames[7:])
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
        self.set_abvega_offset()

        self.make_segment_catalog()
        self.make_null_column()

        self.set_aperture_params()
        self.make_aperture_catalog()
        self.make_extras_catalog()
        self.make_source_catalog()

        # output formatting requested by the JPWG (2020.02.05)
        for colname in self.catalog.colnames:
            if colname in ('xcentroid', 'ycentroid') or 'CI_' in colname:
                self.catalog[colname].info.format = '.4f'
            if 'flux' in colname:
                self.catalog[colname].info.format = '.6e'
            if 'abmag' in colname or 'vegamag' in colname or 'nn_' in colname:
                self.catalog[colname].info.format = '.6f'
            if colname in ('semimajor_sigma', 'semiminor_sigma',
                           'ellipticity', 'orientation', 'sky_orientation'):
                self.catalog[colname].info.format = '.6f'

        return self.catalog
