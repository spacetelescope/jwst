"""
Module to calculate the source catalog.
"""

from collections import OrderedDict
import logging
import warnings

from astropy import __version__ as astropy_version
from astropy.convolution import Gaussian2DKernel
from astropy.nddata.utils import extract_array
from astropy.stats import gaussian_fwhm_to_sigma, SigmaClip
from astropy.table import QTable
import astropy.units as u
from astropy.utils import lazyproperty
import numpy as np
from scipy import __version__ as scipy_version
from scipy import ndimage
from scipy.spatial import cKDTree

from photutils import __version__ as photutils_version
from photutils import Background2D, MedianBackground
from photutils import detect_sources, deblend_sources, source_properties
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from photutils.detection.findstars import _StarFinderKernel
from photutils.utils._wcs_helpers import _pixel_scale_angle_at_skycoord

from jwst import __version__ as jwst_version

from .. import datamodels
from ..datamodels import ImageModel, ABVegaOffsetModel

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class ReferenceData:
    """
    Class for APCORR and ABVEGAOFFSET reference file data needed by
    `SourceCatalogStep`.

    Parameters
    ----------
    model : `ImageModel`
        An `ImageModel` of drizzled image.

    aperture_ee : tuple of 3 int
        The aperture encircled energies to be used for aperture
        photometry.  The values must be 3 strictly-increasing integers.
        Valid values are defined in the APCORR reference files (20, 30,
        40, 50, 60, 70, or 80).

    apcorr_filename : str
        The full path filename of the APCORR reference file.

    abvegaoffset_filename : str
        The full path filename of the ABVEGAOFFSET reference file.

    Attributes
    ----------
    aperture_params : `dict`
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).

    abvega_offset : float
        Offset to convert from AB to Vega magnitudes.  The value
        represents m_AB - m_Vega.
    """

    def __init__(self, model, aperture_ee=(30, 50, 70),
                 apcorr_filename=None, abvegaoffset_filename=None):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model

        self.aperture_ee = self._validate_aperture_ee(aperture_ee)
        self.apcorr_filename = apcorr_filename
        self.abvegaoffset_filename = abvegaoffset_filename

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
        """
        Validate the input ``aperture_ee``.
        """
        aperture_ee = np.array(aperture_ee).astype(int)
        if not np.all(aperture_ee[1:] > aperture_ee[:-1]):
            raise ValueError('aperture_ee values must be strictly '
                             'increasing')
        if len(aperture_ee) != 3:
            raise ValueError('aperture_ee must contain only 3 values')
        if np.any(np.logical_or(aperture_ee <= 0, aperture_ee >= 100)):
            raise ValueError('aperture_ee values must be between 0 and 100')
        return aperture_ee

    @lazyproperty
    def _aperture_ee_table(self):
        """
        Get the encircled energy table for the given instrument
        configuration.
        """
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

        if len(ee_table) == 0:
            raise RuntimeError('APCORR reference file data is missing for '
                               f'{selector}.')

        return ee_table

    def _get_ee_table_row(self, aperture_ee):
        """
        Get the encircled energy row for the input ``aperture_ee``.
        """
        ee_percent = np.round(self._aperture_ee_table['eefraction'] * 100)
        row_mask = (ee_percent == aperture_ee)
        ee_row = self._aperture_ee_table[row_mask]
        if len(ee_row) == 0:
            raise RuntimeError('Aperture encircled energy value of {0} '
                               'appears to be invalid. No matching row '
                               'was found in the APCORR reference file '
                               '{1}'.format(aperture_ee,
                                            self.apcorr_filename))
        if len(ee_row) > 1:
            raise RuntimeError('More than one matching row was found in '
                               'the APCORR reference file {0}'
                               .format(self.apcorr_filename))
        return ee_row

    @lazyproperty
    def aperture_params(self):
        """
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).
        """
        if self.apcorr_filename is None:
            log.warning('APCorrModel reference file was not input. Using '
                        'fallback aperture sizes without any aperture '
                        'corrections.')

            params = {'aperture_radii': np.array((1.0, 2.0, 3.0)),
                      'aperture_corrections': np.array((1.0, 1.0, 1.0)),
                      'aperture_ee': np.array((1, 2, 3)),
                      'bkg_aperture_inner_radius': 5.0,
                      'bkg_aperture_outer_radius': 10.0}
            return params

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
        params['aperture_corrections'] = np.array(apcorrs)

        skyins = np.unique(skyins)
        skyouts = np.unique(skyouts)
        if len(skyins) != 1 or len(skyouts) != 1:
            raise RuntimeError('Expected to find only one value for skyin '
                               'and skyout in the APCORR reference file for '
                               'a given selector.')
        params['bkg_aperture_inner_radius'] = skyins[0]
        params['bkg_aperture_outer_radius'] = skyouts[0]

        return params

    @lazyproperty
    def abvega_offset(self):
        """
        Offset to convert from AB to Vega magnitudes.

        The value represents m_AB - m_Vega.
        """
        if self.abvegaoffset_filename is None:
            log.warning('ABVEGAOFFSET reference file was not input. '
                        'Catalog Vega magnitudes are not correct.')
            return 0.0

        if self.instrument == 'NIRCAM' or self.instrument == 'NIRISS':
            selector = {'filter': self.filtername, 'pupil': self.pupil}
        elif self.instrument == 'MIRI':
            selector = {'filter': self.filtername}
        elif self.instrument == 'FGS':
            selector = {'detector': self.detector}
        else:
            raise RuntimeError(f'{self.instrument} is not a valid instrument')

        abvegaoffset_model = ABVegaOffsetModel(self.abvegaoffset_filename)
        offsets_table = abvegaoffset_model.abvega_offset

        try:
            mask_idx = [offsets_table[key] == value
                        for key, value in selector.items()]
        except KeyError as badkey:
            raise KeyError('{0} not found in ABVEGAOFFSET reference '
                           'file {1}'.format(badkey,
                                             self.abvegaoffset_filename))

        row = offsets_table[np.logical_and.reduce(mask_idx)]

        if len(row) == 0:
            raise RuntimeError('Did not find matching row in ABVEGAOFFSET '
                               'reference file {0}'
                               .format(self.abvegaoffset_filename))
        if len(row) > 1:
            raise RuntimeError('Found more than one matching row in '
                               'ABVEGAOFFSET reference file {0}'
                               .format(self.abvegaoffset_filename))

        abvega_offset = row['abvega_offset'][0]
        log.info('AB to Vega magnitude offset {:.5f}'.format(abvega_offset))
        abvegaoffset_model.close()
        return abvega_offset


class Background:
    """
    Class to estimate a 2D background and background RMS noise in an
    image.

    Parameters
    ----------
    data : 2D `~numpy.ndarray`
        The input 2D array.

    box_size : int or array_like (int)
        The box size along each axis.  If ``box_size`` is a scalar then
        a square box of size ``box_size`` will be used.  If ``box_size``
        has two elements, they should be in ``(ny, nx)`` order.

    mask : array_like (bool), optional
        A boolean mask, with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked data are excluded from calculations.

    Attributes
    ----------
    background : 2D `~numpy.ndimage`
        The estimated 2D background image.

    background_rms : 2D `~numpy.ndimage`
        The estimated 2D background RMS image.
    """

    def __init__(self, data, box_size=100, mask=None):
        self.data = data
        self.box_size = np.asarray(box_size).astype(int)  # must be integer
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

        try:
            bkg = Background2D(self.data, self.box_size,
                               filter_size=filter_size, mask=self.mask,
                               sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)
        except ValueError:
            # use the entire unmasked array
            bkg = Background2D(self.data, self.data.shape,
                               filter_size=filter_size, mask=self.mask,
                               sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator,
                               exclude_percentile=100.)
            log.info('Background could not be estimated in meshes. '
                     'Using the entire unmasked array for background '
                     f'estimation:  bkg_boxsize={self.data.shape}.')

        # apply the coverage mask
        bkg.background *= np.logical_not(self.mask)
        bkg.background_rms *= np.logical_not(self.mask)

        return bkg

    @lazyproperty
    def background(self):
        """
        The 2D background image.
        """
        return self._background2d.background

    @lazyproperty
    def background_rms(self):
        """
        The 2D background RMS image.
        """
        return self._background2d.background_rms


def make_kernel(kernel_fwhm):
    """
    Make a 2D Gaussian smoothing kernel that is used to filter the image
    before thresholding.

    Filtering the image will smooth the noise and maximize detectability
    of objects with a shape similar to the kernel.

    The kernel must have odd sizes in both X and Y, be centered in the
    central pixel, and normalized to sum to 1.

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
    data : 2D `~numpy.ndarray`
        The input 2D array.

    threshold : float
        The data value or pixel-wise data values to be used for the
        detection threshold. A 2D threshold must have the same shape as
        ``data``.

    npixels : int
        The number of connected pixels, each greater than ``threshold``
        that an object must have to be detected.  ``npixels`` must be a
        positive integer.

    kernel : `astropy.convolution.Kernel2D`
        The filtering kernel.  Filtering the image will smooth the noise
        and maximize detectability of objects with a shape similar to
        the kernel.

    mask : array_like of bool, optional
        A boolean mask, with the same shape as the input ``data``, where
        `True` values indicate masked pixels.  Masked pixels will not be
        included in any source.

    deblend : bool, optional
        Whether to deblend overlapping sources.  Source deblending
        requires scikit-image.

    Returns
    -------
    segment_image : `~photutils.segmentation.SegmentationImage` or `None`
        A 2D segmentation image, with the same shape as the input data,
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
    Calculate the total error array.

    The total error includes both background-only error and the source
    Poisson noise.

    Parameters
    ----------
    model : `ImageModel`
        The input `ImageModel`.

    Returns
    -------
    total_error : `~numpy.ndarray`
        The total error array.
    """
    # TODO: Until errors are produced for the level 3 drizzle
    # products, the JPWG has decided that the errors should not be
    # populated.
    return np.zeros_like(model.data)


class SourceCatalog:
    """
    Class for the source catalog.

    Parameters
    ----------
    model : `ImageModel`
        The input `ImageModel`.  The data is assumed to be
        background subtracted.

    segment_image : `~photutils.segmentation.SegmentationImage`
        A 2D segmentation image, with the same shape as the input data,
        where sources are marked by different positive integer values.
        A value of zero is reserved for the background.

    ci_star_thresholds : array-like of 2 floats
        The concentration index thresholds for determining whether
        a source is a star. The first threshold corresponds to the
        concentration index calculated from the smallest and middle
        aperture radii (see ``aperture_params``). The second threshold
        corresponds to the concentration index calculated from the
        middle and largest aperture radii. An object is considered
        extended if both concentration indices are greater than the
        corresponding thresholds, otherwise it is considered a star.

    error : array_like or `~astropy.units.Quantity`, optional
        The total error array corresponding to the input ``data`` array.
        ``error`` is assumed to include *all* sources of error,
        including the Poisson error of the sources (see
        `~photutils.utils.calc_total_error`) .  ``error`` must have the
        same shape as the input ``data``.

    kernel : array-like (2D) or `~astropy.convolution.Kernel2D`, optional
        The 2D array of the kernel used to filter the data prior to
        calculating the source centroid and morphological parameters.
        The kernel should be the same one used in defining the source
        segments.

    kernel_fwhm : float
        The full-width at half-maximum (FWHM) of the 2D Gaussian kernel.

    aperture_params : `dict`
        A dictionary containing the aperture parameters (radii, aperture
        corrections, and background annulus inner and outer radii).

    abvega_offset : float
        Offset to convert from AB to Vega magnitudes.  The value
        represents m_AB - m_Vega.
    """

    def __init__(self, model, segment_img, ci_star_thresholds, error=None,
                 kernel=None, kernel_fwhm=None, aperture_params=None,
                 abvega_offset=0.0):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model  # background was subtracted in SourceDetection

        self.segment_img = segment_img
        if len(ci_star_thresholds) != 2:
            raise ValueError('ci_star_thresholds must contain only 2 '
                             'items')
        self.ci_star_thresholds = ci_star_thresholds
        self.error = error  # total error array
        self.kernel = kernel
        self.kernel_fwhm = kernel_fwhm
        self.aperture_params = aperture_params
        self.abvega_offset = abvega_offset

        self.aperture_ee = aperture_params['aperture_ee']
        self.n_aper = len(self.aperture_ee)
        self.column_desc = {}

        self.wcs = self.model.meta.wcs
        self._xpeak = None
        self._ypeak = None

    def convert_to_jy(self):
        """
        Convert the data and errors from MJy/sr to Jy and convert to
        `~astropy.unit.Quantity` objects.
        """
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
        Convert flux (and error) to AB magnitude (and error).

        Parameters
        ----------
        flux, flux_err : `~astropy.unit.Quantity`
            The input flux and error arrays in units of Jy.

        Returns
        -------
        abmag, abmag_err : `~astropy.ndarray`
            The output AB magnitude and error arrays.
        """
        # ignore RunTimeWarning if flux or flux_err contains NaNs
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)

            abmag = -2.5 * np.log10(flux.value) + 8.9
            abmag_err = 2.5 * np.log10(1.0 + (flux_err.value / flux.value))

            # handle negative fluxes
            idx = flux.value < 0
            abmag[idx] = np.nan
            abmag_err[idx] = np.nan

        return abmag, abmag_err

    @lazyproperty
    def segment_colnames(self):
        """
        A dictionary of the output table column names and descriptions
        for the segment catalog.
        """
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
        Calculate the segment-based source photometry and morphologies.

        The values are set as dynamic attributes.
        """
        source_props = source_properties(self.model.data.astype(float),
                                         self.segment_img,
                                         error=self.error,
                                         filter_kernel=self.kernel,
                                         wcs=self.wcs)

        self._xpeak = source_props.maxval_xpos.value.astype(int)
        self._ypeak = source_props.maxval_ypos.value.astype(int)

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
        """
        An array containing only NaNs.
        """
        values = np.empty(len(self.id))
        values.fill(np.nan)
        return values

    @lazyproperty
    def xypos(self):
        """
        The (x, y) source positions, defined from the segmentation
        image.
        """
        return np.transpose((self.xcentroid, self.ycentroid))

    @lazyproperty
    def _isophotal_abmag(self):
        """
        The isophotal AB magnitude and error.
        """
        return self.convert_flux_to_abmag(self.isophotal_flux,
                                          self.isophotal_flux_err)

    @lazyproperty
    def isophotal_abmag(self):
        """
        The isophotal AB magnitude.
        """
        return self._isophotal_abmag[0]

    @lazyproperty
    def isophotal_abmag_err(self):
        """
        The isophotal AB magnitude error.
        """
        return self._isophotal_abmag[1]

    @lazyproperty
    def isophotal_vegamag(self):
        """
        The isophotal Vega magnitude.
        """
        return self.isophotal_abmag - self.abvega_offset

    @lazyproperty
    def isophotal_vegamag_err(self):
        """
        The isophotal Vega magnitude error.
        """
        return self.isophotal_abmag_err

    @lazyproperty
    def sky_orientation(self):
        """
        The orientation of the source major axis as the position angle
        in degrees measured East of North.
        """
        # NOTE: crpix1 and crpix2 are 1-based values
        skycoord = self.wcs.pixel_to_world(self.model.meta.wcsinfo.crpix1 - 1,
                                           self.model.meta.wcsinfo.crpix2 - 1)
        _, angle = _pixel_scale_angle_at_skycoord(skycoord, self.wcs)

        return (180.0 * u.deg) - angle + self.orientation

    def _make_aperture_colnames(self, name):
        """
        Make the aperture column names.

        There are separate columns for the flux/magnitude for each of
        the encircled energies and the total flux/magnitude.

        Parameters
        ----------
        name : {'flux', 'abmag', 'vegamag'}
            The name type of the column.

        Returns
        -------
        colnames : list of str
            A list of the output column names.
        """
        colnames = []
        for aper_ee in self.aperture_ee:
            basename = f'aper{aper_ee}_{name}'
            colnames.append(basename)
            colnames.append(f'{basename}_err')
        colnames.extend([f'aper_total_{name}', f'aper_total_{name}_err'])

        return colnames

    def _make_aperture_descriptions(self, name):
        """
        Make aperture column descriptions.

        Parameters
        ----------
        name : {'flux', 'abmag', 'vegamag'}
            The name type of the column.

        Returns
        -------
        descriptions : list of str
            A list of the output column descriptions.
        """
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
                    'aperture; calculated only for stars')
        desc.append(f'Total aperture-corrected {ftype2} error based on the '
                    f'{self.aperture_ee[-1]}% encircled energy circular '
                    'aperture; calculated only for stars')

        return desc

    @lazyproperty
    def aperture_flux_colnames(self):
        """
        The aperture flux column names.
        """
        return self._make_aperture_colnames('flux')

    @lazyproperty
    def aperture_flux_descriptions(self):
        """
        The aperture flux column descriptions.
        """
        return self._make_aperture_descriptions('flux')

    @lazyproperty
    def aperture_abmag_colnames(self):
        """
        The aperture AB magnitude column names.
        """
        return self._make_aperture_colnames('abmag')

    @lazyproperty
    def aperture_abmag_descriptions(self):
        """
        The aperture AB magnitude column descriptions.
        """
        return self._make_aperture_descriptions('abmag')

    @lazyproperty
    def aperture_vegamag_colnames(self):
        """
        The aperture Vega magnitude column names.
        """
        return self._make_aperture_colnames('vegamag')

    @lazyproperty
    def aperture_vegamag_descriptions(self):
        """
        The aperture Vega magnitude column descriptions.
        """
        return self._make_aperture_descriptions('vegamag')

    @lazyproperty
    def aperture_colnames(self):
        """
        A dictionary of the output table column names and descriptions
        for the aperture catalog.
        """
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
        Estimate the local background and error using a circular annulus
        aperture.

        The local backround is the sigma-clipped median value in the
        annulus.  The background error is the standard error of the
        median, sqrt(pi / 2N) * std.
        """
        bkg_aper = CircularAnnulus(
            self.xypos, self.aperture_params['bkg_aperture_inner_radius'],
            self.aperture_params['bkg_aperture_outer_radius'])
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
        """
        The aperture local background flux (per pixel).
        """
        return self._aper_local_background[0]

    @lazyproperty
    def aper_bkg_flux_err(self):
        """
        The aperture local background flux error (per pixel).
        """
        return self._aper_local_background[1]

    def set_aperture_properties(self):
        """
        Calculate the aperture photometry.

        The values are set as dynamic attributes.
        """
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
        """
        A dictionary of the output table column names and descriptions
        for the additional catalog values.
        """
        desc = OrderedDict()
        for idx, colname in enumerate(self.ci_colnames):
            desc[colname] = self.ci_colname_descriptions[idx]

        desc['is_star'] = 'Flag indicating whether the source is a star'
        desc['sharpness'] = 'The DAOFind source sharpness statistic'
        desc['roundness'] = 'The DAOFind source roundness statistic'
        desc['nn_dist'] = 'The distance in pixels to the nearest neighbor'
        desc['nn_abmag'] = ('The AB magnitude of the nearest neighbor.  If '
                            'the object is a star it is the total '
                            'aperture-corrected AB magnitude, otherwise it '
                            'is the isophotal AB magnitude.')

        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def _ci_ee_indices(self):
        """
        The EE indicies for the concentration indices.

        The three concentration indices are the difference in
        AB magnitudes between:

            * the smallest and middle aperture radii/EE;  idx = (0, 1)
            * the middle and largest aperture radii/EE; idx (1, 2)
            * the smallest and largest aperture radii/EE; idx (0, 2)
        """
        # NOTE: the EE values are always in increasing order
        return ((0, 1), (1, 2), (0, 2))

    @lazyproperty
    def ci_colnames(self):
        """
        The column names of the three concentration indices.
        """
        return [f'CI_{self.aperture_ee[i]}_{self.aperture_ee[j]}'
                for (i, j) in self._ci_ee_indices]

    @lazyproperty
    def ci_colname_descriptions(self):
        """
        The concentration indicies column descriptions.
        """
        return ['Concentration index calculated as '
                f'{self.aperture_abmag_colnames[2*i]} - '
                f'{self.aperture_abmag_colnames[2*j]}' for (i, j) in
                self._ci_ee_indices]

    @lazyproperty
    def concentration_indices(self):
        """
        A list of concentration indices, calculated as the difference of
        AB magnitudes between:

            * the smallest and middle aperture radii/EE
            * the middle and largest aperture radii/EE
            * the smallest and largest aperture radii/EE
        """
        abmags = [(self.aperture_abmag_colnames[2*i],
                   self.aperture_abmag_colnames[2*j]) for (i, j) in
                  self._ci_ee_indices]
        return [getattr(self, mag1) - getattr(self, mag2)
                for mag1, mag2 in abmags]

    def set_ci_properties(self):
        """
        Set the concentration indices as dynamic attributes.
        """
        for name, value in zip(self.ci_colnames, self.concentration_indices):
            setattr(self, name, value)
        return

    @lazyproperty
    def is_star(self):
        """
        Boolean indicating whether the source is a star.
        """
        mask1 = self.concentration_indices[0] > self.ci_star_thresholds[0]
        mask2 = self.concentration_indices[1] > self.ci_star_thresholds[1]
        return np.logical_not(np.logical_and(mask1, mask2))

    @lazyproperty
    def _daofind_kernel(self):
        """
        The DAOFind kernel.
        """
        kernel = _StarFinderKernel(self.kernel_fwhm, ratio=1.0, theta=0.0,
                                   sigma_radius=1.5, normalize_zerosum=True)
        return kernel

    @lazyproperty
    def _kernel_xcenter(self):
        """
        The DAOFind kernel x center.
        """
        return (self._daofind_kernel.data.shape[1] - 1) // 2

    @lazyproperty
    def _kernel_ycenter(self):
        """
        The DAOFind kernel y center.
        """
        return (self._daofind_kernel.data.shape[0] - 1) // 2

    @lazyproperty
    def _daofind_convolved_data(self):
        """
        The DAOFind convolved data.
        """
        return ndimage.convolve(self.model.data.value,
                                self._daofind_kernel.data, mode='constant',
                                cval=0.0)

    @lazyproperty
    def _daofind_cutout(self):
        """
        3D array containing 2D cutouts centered on each source from the
        input data.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        cutout = []
        for xpeak, ypeak in zip(self._xpeak, self._ypeak):
            cutout.append(extract_array(self.model.data,
                                        self._daofind_kernel.data.shape,
                                        (ypeak, xpeak),
                                        fill_value=0.0))
        return np.array(cutout)  # all cutouts are the same size

    @lazyproperty
    def _daofind_cutout_conv(self):
        """
        3D array containing 2D cutouts centered on each source from the
        DAOFind convolved data.

        The cutout size always matches the size of the DAOFind kernel,
        which has odd dimensions.
        """
        cutout = []
        for xpeak, ypeak in zip(self._xpeak, self._ypeak):
            cutout.append(extract_array(self._daofind_convolved_data,
                                        self._daofind_kernel.data.shape,
                                        (ypeak, xpeak),
                                        fill_value=0.0))
        return np.array(cutout)  # all cutouts are the same size

    @lazyproperty
    def sharpness(self):
        """
        The DAOFind source sharpness statistic.

        The sharpness statistic measures the ratio of the difference
        between the height of the central pixel and the mean of the
        surrounding non-bad pixels to the height of the best fitting
        Gaussian function at that point.

        Stars generally have a ``sharpness`` between 0.2 and 1.0.
        """
        npixels = self._daofind_kernel.npixels - 1  # exclude the peak pixel
        data_masked = self._daofind_cutout * self._daofind_kernel.mask
        data_peak = self._daofind_cutout[:, self._kernel_ycenter,
                                         self._kernel_xcenter]
        conv_peak = self._daofind_cutout_conv[:, self._kernel_ycenter,
                                              self._kernel_xcenter]

        data_mean = ((np.sum(data_masked, axis=(1, 2)) -
                      data_peak) / npixels)

        return (data_peak - data_mean) / conv_peak

    @lazyproperty
    def roundness(self):
        """
        The DAOFind source roundness statistic based on symmetry.

        The roundness characteristic computes the ratio of a measure of
        the bilateral symmetry of the object to a measure of the
        four-fold symmetry of the object.

        "Round" objects have a ``roundness`` close to 0, generally
        between -1 and 1.
        """
        # set the central (peak) pixel to zero
        cutout = self._daofind_cutout_conv.copy()
        cutout[:, self._kernel_ycenter, self._kernel_xcenter] = 0.0

        # calculate the four roundness quadrants
        quad1 = cutout[:, 0:self._kernel_ycenter + 1,
                       self._kernel_xcenter + 1:]
        quad2 = cutout[:, 0:self._kernel_ycenter, 0:self._kernel_xcenter + 1]
        quad3 = cutout[:, self._kernel_ycenter:, 0:self._kernel_xcenter]
        quad4 = cutout[:, self._kernel_ycenter + 1:, self._kernel_xcenter:]

        axis = (1, 2)
        sum2 = (-quad1.sum(axis=axis) + quad2.sum(axis=axis) -
                quad3.sum(axis=axis) + quad4.sum(axis=axis))
        sum2[sum2 == 0] = 0.0

        sum4 = np.abs(cutout).sum(axis=axis)
        sum4[sum4 == 0] = np.nan

        return 2.0 * sum2 / sum4

    @lazyproperty
    def _ckdtree_query(self):
        """
        The distance in pixels to the nearest neighbor and its index.
        """
        if self.xypos.shape[0] == 1:  # only one detected source
            return [np.nan], [np.nan]
        else:
            tree = cKDTree(self.xypos)
            qdist, qidx = tree.query(self.xypos, k=2)
            return np.transpose(qdist)[1], np.transpose(qidx)[1]

    @lazyproperty
    def nn_dist(self):
        """
        The distance in pixels to the nearest neighbor.
        """
        return self._ckdtree_query[0] * u.pixel

    @lazyproperty
    def nn_abmag(self):
        """
        The AB magnitude of the nearest neighbor.

        If the nearest-neighbor is a star, then its aperture-corrected
        total AB magnitude is used.  Otherwise, its isophotal AB
        magnitude is used.
        """
        if len(self._ckdtree_query[1]) == 1:  # only one detected source
            return np.nan

        nn_abmag = self.aper_total_abmag[self._ckdtree_query[1]]
        nn_isomag = self.isophotal_abmag[self._ckdtree_query[1]]
        return np.where(np.isnan(nn_abmag), nn_isomag, nn_abmag)

    @lazyproperty
    def _not_star(self):
        """
        Whether the source is not a star.
        """
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
        The aperture-corrected total flux for sources that are stars,
        based on the flux in largest aperture.
        """
        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux = (self.aperture_params['aperture_corrections'][idx] *
                getattr(self, self.aperture_flux_colnames[idx*2]))
        flux[self._not_star] = np.nan

        return flux

    @lazyproperty
    def aper_total_flux_err(self):
        """
        The aperture-corrected total flux error for sources that are
        stars, based on the flux in largest aperture.
        """
        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux_err = (self.aperture_params['aperture_corrections'][idx] *
                    getattr(self, self.aperture_flux_colnames[idx*2 + 1]))
        flux_err[self._not_star] = np.nan

        return flux_err

    @lazyproperty
    def _abmag_total(self):
        """
        The total AB magnitude and error.
        """
        return self.convert_flux_to_abmag(self.aper_total_flux,
                                          self.aper_total_flux_err)

    @lazyproperty
    def aper_total_abmag(self):
        """
        The total AB magnitude.
        """
        return self._abmag_total[0]

    @lazyproperty
    def aper_total_abmag_err(self):
        """
        The total AB magnitude error.
        """
        return self._abmag_total[1]

    @lazyproperty
    def aper_total_vegamag(self):
        """
        The total Vega magnitude.
        """
        return self.aper_total_abmag - self.abvega_offset

    @lazyproperty
    def aper_total_vegamag_err(self):
        """
        The total Vega magnitude error.
        """
        return self.aper_total_abmag_err

    @lazyproperty
    def colnames(self):
        """
        The column name order for the final source catalog.
        """
        colnames = self.segment_colnames[0:4]
        colnames.extend(self.aperture_colnames)
        colnames.extend(self.extras_colnames)
        colnames.extend(self.segment_colnames[4:])
        return colnames

    @staticmethod
    def format_columns(catalog):
        """
        Format the values in the output catalog.

        Parameters
        ----------
        catalog : `~astropy.table.Table`
            The catalog to format.

        Returns
        -------
        result : `~astropy.table.Table`
            The formated catalog.
        """
        # output formatting requested by the JPWG (2020.02.05)
        for colname in catalog.colnames:
            if colname in ('xcentroid', 'ycentroid') or 'CI_' in colname:
                catalog[colname].info.format = '.4f'
            if 'flux' in colname:
                catalog[colname].info.format = '.6e'
            if ('abmag' in colname or 'vegamag' in colname or 'nn_' in colname
                    or colname in ('sharpness', 'roundness')):
                catalog[colname].info.format = '.6f'
            if colname in ('semimajor_sigma', 'semiminor_sigma',
                           'ellipticity', 'orientation', 'sky_orientation'):
                catalog[colname].info.format = '.6f'

        return catalog

    @lazyproperty
    def catalog_metadata(self):
        """
        The catalog metadata, include package version numbers.
        """
        meta = {}
        meta['jwst version'] = jwst_version
        meta['numpy version'] = np.__version__
        meta['scipy version'] = scipy_version
        meta['astropy version'] = astropy_version
        meta['photutils version'] = photutils_version
        meta['aperture_params'] = self.aperture_params
        meta['abvega_offset'] = self.abvega_offset
        return meta

    @lazyproperty
    def catalog(self):
        """
        The final source catalog.
        """
        self.convert_to_jy()
        self.set_segment_properties()
        self.set_aperture_properties()
        self.set_ci_properties()

        catalog = QTable()
        for column in self.colnames:
            catalog[column] = getattr(self, column)
            catalog[column].info.description = self.column_desc[column]

        # TODO: Until errors are produced for the level 3 drizzle
        # products, the JPWG has decided that the errors should not be
        # populated.
        for colname in catalog.colnames:
            if colname.endswith('_err'):
                catalog[colname][:] = self.null_column

        catalog = self.format_columns(catalog)
        catalog.meta.update(self.catalog_metadata)

        return catalog
