from collections import OrderedDict
import logging

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip
from astropy.table import QTable
import astropy.units as u
from astropy.utils import lazyproperty
import numpy as np
from scipy.spatial import cKDTree

from photutils import Background2D, MedianBackground
from photutils import detect_sources, deblend_sources, source_properties
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.utils._wcs_helpers import _pixel_scale_angle_at_skycoord

from .. import datamodels
from ..datamodels import ImageModel, ABVegaOffsetModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ReferenceData:
    def __init__(self, model, aperture_ee=(30, 50, 70),
                 apcorr_filename=None, abvega_offset_filename=None):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model

        self.aperture_ee = self._validate_aperture_ee(aperture_ee)
        if apcorr_filename is None:
            raise ValueError('apcorr_filename must be input')
        self.apcorr_filename = apcorr_filename
        self.abvega_offset_filename = abvega_offset_filename

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

    @staticmethod
    def _validate_aperture_ee(aperture_ee):
        aperture_ee = np.array(aperture_ee).astype(int)

        if not np.all(aperture_ee[1:] > aperture_ee[:-1]):
            raise ValueError('aperture_ee values must be strictly '
                             'increasing')
        if len(aperture_ee) != 3:
            raise ValueError('aperture_ee must contain only 3 values')
        if np.any(np.logical_or(aperture_ee <= 0, aperture_ee >= 100)):
            raise ValueError('aperture_ee values must be between 0 and 100')
        if np.any(aperture_ee[1:] < aperture_ee[:-1]):
            raise ValueError('aperture_ee values must be in increasing order')

        return aperture_ee

    @lazyproperty
    def _aperture_ee_table(self):
        if self.instrument == 'NIRCAM' or self.instrument == 'NIRISS':
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername, 'subarray': self.subarray}
        elif self.instrument == 'FGS':
            selector = None
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        apcorr_model = datamodels.open(self.apcorr_filename)
        apcorr = apcorr_model.apcorr_table
        if selector is None:  # FGS
            ee_table = apcorr
        else:
            mask_idx = [apcorr[key] == value
                        for key, value in selector.items()]
            ee_table = apcorr[np.logical_and.reduce(mask_idx)]

        return ee_table

    def _get_ee_table_row(self, aperture_ee):
        ee_percent = np.round(self._aperture_ee_table['eefraction'] * 100)
        row_mask = (ee_percent == aperture_ee)
        ee_row = self._aperture_ee_table[row_mask]
        if len(ee_row) == 0:
            raise RuntimeError('Aperture encircled energy value of {0} '
                               'appears to be invalid. No matching row '
                               'was found in the apcorr reference file '
                               '{1}'.format(aperture_ee,
                                            self.apcorr_filename))
        if len(ee_row) > 1:
            raise RuntimeError('More than one matching row was found in '
                               'the apcorr reference file {0}'
                               .format(self.apcorr_filename))
        return ee_row

    @lazyproperty
    def aperture_params(self):
        params = {}

        radii = []
        apcorrs = []
        skyins = []
        skyouts = []
        for aper_ee in self.aperture_ee:
            row = self._get_ee_table_row(aper_ee)
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
        params['aperture_ee'] = self.aperture_ee
        params['aperture_radii'] = np.array(radii)
        params['aperture_corrs'] = np.array(apcorrs)

        skyins = np.unique(skyins)
        skyouts = np.unique(skyouts)
        if len(skyins) != 1 or len(skyouts) != 1:
            raise RuntimeError('Expected to find only one value for skyin '
                               'and skyout in the apcorr reference file for '
                               'a given selector.')
        params['bkg_aperture_inner'] = skyins[0]
        params['bkg_aperture_outer'] = skyouts[0]

        return params

    @lazyproperty
    def abvega_offset(self):
        if self.instrument == 'NIRCAM' or self.instrument == 'NIRISS':
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername}
        elif self.instrument == 'FGS':
            selector = {'detector': self.detector}
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        if self.abvega_offset_filename is None:
            log.info('ABVegaOffsetModel reference file was not input -- '
                     'catalog Vega magnitudes are not correct.')
            return 0.0

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

        abvega_offset = row['abvega_offset'][0]
        log.info('AB to Vega magnitude offset {:.5f}'.format(abvega_offset))

        return abvega_offset


class Background:
    def __init__(self, data, bkg_boxsize=100, mask=None):
        self.data = data
        self.bkg_boxsize = np.int(bkg_boxsize)  # must be integer
        self.mask = mask

    @lazyproperty
    def _background2d(self):
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
        bkg = Background2D(self.data, self.bkg_boxsize,
                           filter_size=filter_size, mask=self.mask,
                           sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        # apply the coverage mask
        bkg.background *= np.logical_not(self.mask)
        bkg.background_rms *= np.logical_not(self.mask)

        return bkg

    @lazyproperty
    def background(self):
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        return self._background2d.background_rms


def make_kernel(kernel_fwhm):
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

    sigma = kernel_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma)
    kernel.normalize(mode='integral')

    return kernel


def make_segment_img(data, threshold, npixels=5.0, kernel=None, mask=None,
                     deblend=False):
    """
    Detect sources in an image, including deblending.

    Parameters
    ----------
    data : `ImageModel`
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
    segm = detect_sources(data, threshold, npixels, filter_kernel=kernel,
                          mask=mask, connectivity=connectivity)

    # segm=None for photutils >= 0.7
    # segm.nlabels=0 for photutils < 0.7
    if segm is None or segm.nlabels == 0:
        return None

    # source deblending requires scikit-image
    if deblend:
        nlevels = 32
        contrast = 0.001
        mode = 'exponential'
        segm = deblend_sources(data, segm, npixels=npixels,
                               filter_kernel=kernel, nlevels=nlevels,
                               contrast=contrast, mode=mode,
                               connectivity=connectivity, relabel=True)

    return segm


def calc_total_error(model):
    """
    Calculate the total error (background plus source Poisson error).
    """

    # TODO: Until errors are produced for the level 3 drizzle
    # products, the JPWG has decided that the errors should not be
    # populated.
    return np.zeros_like(model.data)


class SourceCatalog:
    def __init__(self, model, segment_img, error=None, kernel=None,
                 aperture_params=None, abvega_offset=0.0):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model  # background was subtracted in SourceDetection

        self.segment_img = segment_img
        self.error = error  # total error array
        self.kernel = kernel
        self.aperture_params = aperture_params
        self.abvega_offset = abvega_offset

        self.aperture_ee = aperture_params['aperture_ee']
        self.n_aper = len(self.aperture_ee)
        self.column_desc = {}

        # self.wcs = self.model.meta.wcs  # gWCS
        self.wcs = self.model.get_fits_wcs()  # FITS WCS

    def convert_to_jy(self):
        """ convert from MJy/sr to Jy"""

        if self.model.meta.bunit_data != 'MJy/sr':
            raise ValueError('data is expected to be in units of MJy/sr')
        self.model.data *= (1.e6 *
                            self.model.meta.photometry.pixelarea_steradians)
        self.model.meta.bunit_data = 'Jy'

        unit = u.Jy
        self.model.data <<= unit
        self.error <<= unit

        return

    @staticmethod
    def convert_flux_to_abmag(flux, flux_err):
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

    @lazyproperty
    def segment_colnames(self):
        # define segment column names and descriptions
        desc = OrderedDict()
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

        self.column_desc.update(desc)

        return list(desc.keys())

    def set_segment_properties(self):
        """
        Create a catalog of source photometry and morphologies based on
        segmentation.
        """

        source_props = source_properties(self.model.data.astype(float),
                                         self.segment_img,
                                         error=self.error,
                                         filter_kernel=self.kernel,
                                         wcs=self.wcs)

        # rename some columns in the output catalog
        prop_names = {}
        prop_names['isophotal_flux'] = 'source_sum'
        prop_names['isophotal_flux_err'] = 'source_sum_err'
        prop_names['isophotal_area'] = 'area'
        prop_names['semimajor_sigma'] = 'semimajor_axis_sigma'
        prop_names['semiminor_sigma'] = 'semiminor_axis_sigma'

        for column in self.segment_colnames:
            # define the property name
            prop_name = prop_names.get(column, column)

            try:
                value = getattr(source_props, prop_name)
            except AttributeError:
                value = getattr(self, prop_name)

            setattr(self, column, value)

        return

    @lazyproperty
    def null_column(self):
        values = np.empty(len(self.id))
        values.fill(np.nan)
        return values

    @lazyproperty
    def xypos(self):
        # source positions
        return np.transpose((self.xcentroid, self.ycentroid))

    @lazyproperty
    def _isophotal_abmag(self):
        return self.convert_flux_to_abmag(self.isophotal_flux,
                                          self.isophotal_flux_err)

    @lazyproperty
    def isophotal_abmag(self):
        return self._isophotal_abmag[0]

    @lazyproperty
    def isophotal_abmag_err(self):
        return self._isophotal_abmag[1]

    @lazyproperty
    def isophotal_vegamag(self):
        return self.isophotal_abmag - self.abvega_offset

    @lazyproperty
    def isophotal_vegamag_err(self):
        return self.isophotal_abmag_err

    @lazyproperty
    def sky_orientation(self):
        # Define the orientation position angle (E of N)
        # NOTE: crpix1 and crpix2 are 1-based values
        skycoord = self.wcs.pixel_to_world(self.model.meta.wcsinfo.crpix1 - 1,
                                           self.model.meta.wcsinfo.crpix2 - 1)
        _, angle = _pixel_scale_angle_at_skycoord(skycoord, self.wcs)

        return (180.0 * u.deg) - angle + self.orientation

    def _make_aperture_colnames(self, name):
        colnames = []
        for aper_ee in self.aperture_ee:
            basename = f'aper{aper_ee}_{name}'
            colnames.append(basename)
            colnames.append(f'{basename}_err')
        colnames.extend([f'aper_total_{name}', f'aper_total_{name}_err'])

        return colnames

    def _make_aperture_descriptions(self, name):
        if name == 'flux':
            ftype = 'Flux'
            ftype2 = 'flux'
        elif name == 'abmag':
            ftype = ftype2 = 'AB magnitude'
        elif name == 'vegamag':
            ftype = ftype2 = 'Vega magnitude'

        desc = []
        for aper_ee in self.aperture_ee:
            desc.append(f'{ftype} within the {aper_ee}% encircled energy '
                        'circular aperture')
            desc.append(f'{ftype} error within the {aper_ee}% encircled '
                        'energy circular aperture')

        desc.append(f'Total aperture-corrected {ftype2} based on the '
                    f'{self.aperture_ee[-1]}% encircled energy circular '
                    'aperture - calculated only for stars')
        desc.append(f'Total aperture-corrected {ftype2} error based on the '
                    f'{self.aperture_ee[-1]}% encircled energy circular '
                    'aperture - calculated only for stars')

        return desc

    @lazyproperty
    def aperture_flux_colnames(self):
        return self._make_aperture_colnames('flux')

    @lazyproperty
    def aperture_flux_descriptions(self):
        return self._make_aperture_descriptions('flux')

    @lazyproperty
    def aperture_abmag_colnames(self):
        return self._make_aperture_colnames('abmag')

    @lazyproperty
    def aperture_abmag_descriptions(self):
        return self._make_aperture_descriptions('abmag')

    @lazyproperty
    def aperture_vegamag_colnames(self):
        return self._make_aperture_colnames('vegamag')

    @lazyproperty
    def aperture_vegamag_descriptions(self):
        return self._make_aperture_descriptions('vegamag')

    @lazyproperty
    def aperture_colnames(self):
        desc = OrderedDict()
        desc['aper_bkg_flux'] = ('The local background value calculated as '
                                 'the sigma-clipped median value in the '
                                 'background annulus aperture')
        desc['aper_bkg_flux_err'] = ('The standard error of the '
                                     'sigma-clipped median background value')

        for idx, colname in enumerate(self.aperture_flux_colnames):
            desc[colname] = self.aperture_flux_descriptions[idx]
        for idx, colname in enumerate(self.aperture_abmag_colnames):
            desc[colname] = self.aperture_abmag_descriptions[idx]
        for idx, colname in enumerate(self.aperture_vegamag_colnames):
            desc[colname] = self.aperture_vegamag_descriptions[idx]

        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def _aper_local_background(self):
        """
        Estimate the local background as the sigma-clipped median value
        in the input aperture.
        """

        bkg_aper = CircularAnnulus(self.xypos,
                                   self.aperture_params['bkg_aperture_inner'],
                                   self.aperture_params['bkg_aperture_outer'])
        bkg_aper_masks = bkg_aper.to_mask(method='center')
        sigclip = SigmaClip(sigma=3)

        nvalues = []
        bkg_median = []
        bkg_std = []
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

    @lazyproperty
    def aper_bkg_flux(self):
        return self._aper_local_background[0]

    @lazyproperty
    def aper_bkg_flux_err(self):
        return self._aper_local_background[1]

    def set_aperture_properties(self):
        apertures = [CircularAperture(self.xypos, radius) for radius in
                     self.aperture_params['aperture_radii']]
        aper_phot = aperture_photometry(self.model.data, apertures,
                                        error=self.error)

        for i, aperture in enumerate(apertures):
            flux_col = f'aperture_sum_{i}'
            flux_err_col = f'aperture_sum_err_{i}'

            # subtract the local background measured in the annulus
            aper_phot[flux_col] -= (self.aper_bkg_flux * aperture.area)

            flux = aper_phot[flux_col]
            flux_err = aper_phot[flux_err_col]
            abmag, abmag_err = self.convert_flux_to_abmag(flux, flux_err)
            vegamag = abmag - self.abvega_offset
            vegamag_err = abmag_err

            idx0 = 2 * i
            idx1 = (2 * i) + 1
            setattr(self, self.aperture_flux_colnames[idx0], flux)
            setattr(self, self.aperture_flux_colnames[idx1], flux_err)
            setattr(self, self.aperture_abmag_colnames[idx0], abmag)
            setattr(self, self.aperture_abmag_colnames[idx1], abmag_err)
            setattr(self, self.aperture_vegamag_colnames[idx0], vegamag)
            setattr(self, self.aperture_vegamag_colnames[idx1], vegamag_err)

    @lazyproperty
    def extras_colnames(self):
        desc = OrderedDict()

        ci_desc = ('Concentration index calculated as '
                   f'{self.aperture_abmag_colnames[0]} - '
                   f'{self.aperture_abmag_colnames[4]}')
        desc['concentration_index'] = ci_desc

        desc['is_star'] = 'Flag indicating whether the source is a star'
        desc['sharpness'] = 'The DAOFind source sharpness statistic'
        desc['roundness'] = 'The DAOFind source roundness statistic'
        desc['nn_dist'] = 'The distance in pixels to the nearest neighbor'
        desc['nn_abmag'] = 'The AB magnitude of the nearest neighbor'

        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def ci_colname(self):
        return f'CI_{self.aperture_ee[0]}_{self.aperture_ee[2]}'

    @lazyproperty
    def concentration_index(self):
        abmag1 = self.aperture_abmag_colnames[0]  # ee_fraction[0]
        abmag2 = self.aperture_abmag_colnames[4]  # ee_fraction[2]
        return getattr(self, abmag1) - getattr(self, abmag2)

    @lazyproperty
    def is_star(self):
        """
        Calculate a flag based on whether the source is a star.

        2020.03.02: The JWST Photometry Working Group has not determined
        the criteria for this flag.
        """

        # TODO: need algorithm for this flag
        #is_star = np.random.randint(2, size=len(self.id))
        #return is_star.astype(bool)
        return self.null_column

    @lazyproperty
    def sharpness(self):
        return self.null_column

    @lazyproperty
    def roundness(self):
        return self.null_column

    @lazyproperty
    def _ckdtree_query(self):
        """
        Calculate the distance in pixels to the nearest neighbor and
        its AB magnitude.

        If the nearest-neighbor is a star, then its aperture-corrected
        total AB magnitude is used.  Otherwise, its isophotal AB
        magnitude is used.
        """

        tree = cKDTree(self.xypos)
        qdist, qidx = tree.query(self.xypos, k=2)
        return np.transpose(qdist)[1], np.transpose(qidx)[1]

    @lazyproperty
    def nn_dist(self):
        return self._ckdtree_query[0] * u.pixel

    @lazyproperty
    def nn_abmag(self):
        nn_abmag = self.aper_total_abmag[self._ckdtree_query[1]]
        nn_isomag = self.isophotal_abmag[self._ckdtree_query[1]]
        return np.where(np.isnan(nn_abmag), nn_isomag, nn_abmag)

    @lazyproperty
    def _not_star(self):
        if np.isnan(self.is_star[0]):
            # TODO: remove this when ``is_star`` values are available.
            # this array is all False:
            is_star = np.zeros(len(self.id), dtype=bool)
        else:
            is_star = self.is_star

        return np.logical_not(is_star)

    @lazyproperty
    def aper_total_flux(self):
        """
        Calculate the aperture-corrected total flux for sources that are
        stars.

        The aperture fluxes in the largest aperture (i.e., the largest
        encircled energy) are corrected to total flux.
        """

        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux = (self.aperture_params['aperture_corrs'][idx] *
                getattr(self, self.aperture_flux_colnames[idx*2]))
        flux[self._not_star] = np.nan

        return flux

    @lazyproperty
    def aper_total_flux_err(self):
        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux_err = (self.aperture_params['aperture_corrs'][idx] *
                    getattr(self, self.aperture_flux_colnames[idx*2 + 1]))
        flux_err[self._not_star] = np.nan

        return flux_err

    @lazyproperty
    def _abmag_total(self):
        return self.convert_flux_to_abmag(self.aper_total_flux,
                                          self.aper_total_flux_err)

    @lazyproperty
    def aper_total_abmag(self):
        return self._abmag_total[0]

    @lazyproperty
    def aper_total_abmag_err(self):
        return self._abmag_total[1]

    @lazyproperty
    def aper_total_vegamag(self):
        return self.aper_total_abmag - self.abvega_offset

    @lazyproperty
    def aper_total_vegamag_err(self):
        return self.aper_total_abmag_err

    @lazyproperty
    def colnames(self):
        # Define the column name order for the final source catalog
        colnames = self.segment_colnames[0:4]
        colnames.extend(self.aperture_colnames)
        colnames.extend(self.extras_colnames)
        colnames.extend(self.segment_colnames[4:])

        return colnames

    @staticmethod
    def format_columns(catalog):
        # output formatting requested by the JPWG (2020.02.05)
        for colname in catalog.colnames:
            if colname in ('xcentroid', 'ycentroid') or 'CI_' in colname:
                catalog[colname].info.format = '.4f'
            if 'flux' in colname:
                catalog[colname].info.format = '.6e'
            if 'abmag' in colname or 'vegamag' in colname or 'nn_' in colname:
                catalog[colname].info.format = '.6f'
            if colname in ('semimajor_sigma', 'semiminor_sigma',
                           'ellipticity', 'orientation', 'sky_orientation'):
                catalog[colname].info.format = '.6f'

        return catalog

    @lazyproperty
    def catalog(self):
        #self.subtract_background()
        self.convert_to_jy()
        self.set_segment_properties()
        self.set_aperture_properties()

        renamed = {}
        renamed['concentration_index'] = self.ci_colname

        catalog = QTable()
        for column in self.colnames:
            colname = renamed.get(column, column)
            catalog[colname] = getattr(self, column)
            catalog[colname].info.description = self.column_desc[column]

        # TODO: Until errors are produced for the level 3 drizzle
        # products, the JPWG has decided that the errors should not be
        # populated.
        for colname in catalog.colnames:
            if colname.endswith('_err'):
                catalog[colname][:] = self.null_column

        catalog = self.format_columns(catalog)

        return catalog
