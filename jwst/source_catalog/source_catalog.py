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
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import (detect_sources, deblend_sources,
                                    SourceCatalog)
from photutils.aperture import (CircularAperture, CircularAnnulus,
                                aperture_photometry)
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

        if self.model.meta.resample.pixel_scale_ratio is not None:
            # pixel_scale_ratio is the ratio of the resampled to the native
            # pixel scale (values < 1 have smaller resampled pixels)
            pixel_scale_ratio = self.model.meta.resample.pixel_scale_ratio
        else:
            log.warning('model.meta.resample.pixel_scale_ratio was not '
                        'found. Assuming the native detector pixel scale '
                        '(i.e., pixel_scale_ratio = 1)')
            pixel_scale_ratio = 1.0

        params['aperture_ee'] = self.aperture_ee
        params['aperture_radii'] = np.array(radii) / pixel_scale_ratio
        params['aperture_corrections'] = np.array(apcorrs)

        skyins = np.unique(skyins)
        skyouts = np.unique(skyouts)
        if len(skyins) != 1 or len(skyouts) != 1:
            raise RuntimeError('Expected to find only one value for skyin '
                               'and skyout in the APCORR reference file for '
                               'a given selector.')
        params['bkg_aperture_inner_radius'] = skyins[0] / pixel_scale_ratio
        params['bkg_aperture_outer_radius'] = skyouts[0] / pixel_scale_ratio

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
    if segm is None:
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


class JWSTSourceCatalog:
    """
    Class for the JWST source catalog.

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

    Notes
    -----
    ``model.err`` is assumed to be the total error array corresponding
    to the input science ``model.data`` array. It is assumed to include
    *all* sources of error, including the Poisson error of the sources,
    and have the same shape and units as the science data array.
    """

    def __init__(self, model, segment_img, ci_star_thresholds, kernel=None,
                 kernel_fwhm=None, aperture_params=None, abvega_offset=0.0):

        if not isinstance(model, ImageModel):
            raise ValueError('The input model must be a ImageModel.')
        self.model = model  # background was subtracted in SourceDetection

        self.segment_img = segment_img
        if len(ci_star_thresholds) != 2:
            raise ValueError('ci_star_thresholds must contain only 2 '
                             'items')
        self.ci_star_thresholds = ci_star_thresholds
        self.kernel = kernel
        self.kernel_sigma = kernel_fwhm * gaussian_fwhm_to_sigma

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
        self.model.err <<= unit

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
        desc['label'] = 'Unique source identification label number'
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
        segm_cat = SourceCatalog(self.model.data, self.segment_img,
                                 error=self.model.err, kernel=self.kernel,
                                 wcs=self.wcs)

        self._xpeak = segm_cat.maxval_xindex
        self._ypeak = segm_cat.maxval_yindex

        # rename some columns in the output catalog
        prop_names = {}
        prop_names['isophotal_flux'] = 'segment_flux'
        prop_names['isophotal_flux_err'] = 'segment_fluxerr'
        prop_names['isophotal_area'] = 'area'

        for column in self.segment_colnames:
            # define the property name
            prop_name = prop_names.get(column, column)
            try:
                value = getattr(segm_cat, prop_name)
            except AttributeError:
                value = getattr(self, prop_name)
            setattr(self, column, value)

        return

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
                    'aperture; should be used only for unresolved sources.')
        desc.append(f'Total aperture-corrected {ftype2} error based on the '
                    f'{self.aperture_ee[-1]}% encircled energy circular '
                    'aperture; should be used only for unresolved sources.')

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
                                        error=self.model.err)

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

        desc['is_extended'] = 'Flag indicating whether the source is extended'
        desc['sharpness'] = 'The DAOFind source sharpness statistic'
        desc['roundness'] = 'The DAOFind source roundness statistic'
        desc['nn_label'] = 'The label number of the nearest neighbor'
        desc['nn_dist'] = 'The distance in pixels to the nearest neighbor'
        self.column_desc.update(desc)

        return list(desc.keys())

    @lazyproperty
    def _ci_ee_indices(self):
        """
        The EE indicies for the concentration indices.
        """
        # NOTE: the EE values are always in increasing order
        return ((0, 1), (1, 2), (0, 2))

    @lazyproperty
    def ci_colnames(self):
        """
        The column names of the three concentration indices.
        """
        return [f'CI_{self.aperture_ee[j]}_{self.aperture_ee[i]}'
                for (i, j) in self._ci_ee_indices]

    @lazyproperty
    def ci_colname_descriptions(self):
        """
        The concentration indicies column descriptions.
        """
        return ['Concentration index calculated as '
                f'({self.aperture_flux_colnames[2*j]} / '
                f'{self.aperture_flux_colnames[2*i]})' for (i, j) in
                self._ci_ee_indices]

    @lazyproperty
    def concentration_indices(self):
        """
        A list of concentration indices, calculated as the flux
        ratios of:

            * the middle / smallest aperture radii/EE,
              e.g., CI_50_30 = aper50_flux / aper30_flux
            * the largest / middle aperture radii/EE,
              e.g., CI_70_50 = aper70_flux / aper50_flux
            * the largest / smallest aperture radii/EE,
              e.g., CI_70_30 = aper70_flux / aper30_flux
        """
        fluxes = [(self.aperture_flux_colnames[2 * j],
                   self.aperture_flux_colnames[2 * i]) for (i, j) in
                  self._ci_ee_indices]
        return [getattr(self, flux1).value / getattr(self, flux2).value
                for flux1, flux2 in fluxes]

    def set_ci_properties(self):
        """
        Set the concentration indices as dynamic attributes.
        """
        for name, value in zip(self.ci_colnames, self.concentration_indices):
            setattr(self, name, value)
        return

    @lazyproperty
    def is_extended(self):
        """
        Boolean indicating whether the source is extended.
        """
        mask1 = self.concentration_indices[0] > self.ci_star_thresholds[0]
        mask2 = self.concentration_indices[1] > self.ci_star_thresholds[1]
        return np.logical_and(mask1, mask2)

    @lazyproperty
    def _kernel_size(self):
        """
        The DAOFind kernel size (in both x and y dimensions).
        """
        # always odd
        return 2 * int(max(2.0, 1.5 * self.kernel_sigma)) + 1

    @lazyproperty
    def _kernel_center(self):
        """
        The DAOFind kernel x/y center.
        """
        return (self._kernel_size - 1) // 2

    @lazyproperty
    def _kernel_mask(self):
        """
        The DAOFind kernel circular mask.
        NOTE: 1=good pixels, 0=masked pixels
        """
        yy, xx = np.mgrid[0:self._kernel_size, 0:self._kernel_size]
        radius = np.sqrt((xx - self._kernel_center) ** 2
                         + (yy - self._kernel_center) ** 2)
        return (radius <= max(2.0, 1.5 * self.kernel_sigma)).astype(int)

    @lazyproperty
    def _daofind_kernel(self):
        """
        The DAOFind kernel, a 2D circular Gaussian normalized to have
        zero sum.
        """
        size = self._kernel_size
        kernel = Gaussian2DKernel(self.kernel_sigma, x_size=size,
                                  y_size=size).array
        kernel /= np.max(kernel)
        kernel *= self._kernel_mask

        # normalize the kernel to zero sum
        npixels = self._kernel_mask.sum()
        denom = np.sum(kernel**2) - (np.sum(kernel)**2 / npixels)
        return (((kernel - (kernel.sum() / npixels)) / denom)
                * self._kernel_mask)

    @lazyproperty
    def _daofind_convolved_data(self):
        """
        The DAOFind convolved data.
        """
        return ndimage.convolve(self.model.data.value,
                                self._daofind_kernel, mode='constant',
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
                                        self._daofind_kernel.shape,
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
                                        self._daofind_kernel.shape,
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
        npixels = self._kernel_mask.sum() - 1  # exclude the peak pixel
        data_masked = self._daofind_cutout * self._kernel_mask
        data_peak = self._daofind_cutout[:, self._kernel_center,
                                         self._kernel_center]
        conv_peak = self._daofind_cutout_conv[:, self._kernel_center,
                                              self._kernel_center]

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
        cutout[:, self._kernel_center, self._kernel_center] = 0.0

        # calculate the four roundness quadrants
        quad1 = cutout[:, 0:self._kernel_center + 1, self._kernel_center + 1:]
        quad2 = cutout[:, 0:self._kernel_center, 0:self._kernel_center + 1]
        quad3 = cutout[:, self._kernel_center:, 0:self._kernel_center]
        quad4 = cutout[:, self._kernel_center + 1:, self._kernel_center:]

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
    def nn_label(self):
        """
        The label number of the nearest neighbor.
        """
        if len(self._ckdtree_query[1]) == 1:  # only one detected source
            return np.nan
        return self.label[self._ckdtree_query[1]]

    @lazyproperty
    def nn_dist(self):
        """
        The distance in pixels to the nearest neighbor.
        """
        # self._ckdtree_query[0] is NaN if only one detected source
        return self._ckdtree_query[0] * u.pixel

    @lazyproperty
    def aper_total_flux(self):
        """
        The aperture-corrected total flux for sources, based on the flux
        in largest aperture.

        The aperture-corrected total flux should be used only for
        unresolved sources.
        """
        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux = (self.aperture_params['aperture_corrections'][idx] *
                getattr(self, self.aperture_flux_colnames[idx * 2]))
        return flux

    @lazyproperty
    def aper_total_flux_err(self):
        """
        The aperture-corrected total flux error for sources,
        based on the flux in largest aperture.

        The aperture-corrected total flux error should be used only for
        unresolved sources.
        """
        idx = self.n_aper - 1  # apcorr for the largest EE (largest radius)
        flux_err = (self.aperture_params['aperture_corrections'][idx] *
                    getattr(self, self.aperture_flux_colnames[idx * 2 + 1]))
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
        The aperture-corrected total AB magnitude.

        The aperture-corrected total magnitude should be used only for
        unresolved sources.
        """
        return self._abmag_total[0]

    @lazyproperty
    def aper_total_abmag_err(self):
        """
        The aperture-corrected total AB magnitude error.

        The aperture-corrected total magnitude error should be used only
        for unresolved sources.
        """
        return self._abmag_total[1]

    @lazyproperty
    def aper_total_vegamag(self):
        """
        The aperture-corrected total Vega magnitude.

        The aperture-corrected total magnitude should be used only for
        unresolved sources.
        """
        return self.aper_total_abmag - self.abvega_offset

    @lazyproperty
    def aper_total_vegamag_err(self):
        """
        The aperture-corrected total Vega magnitude error.

        The aperture-corrected total magnitude error should be used only
        for unresolved sources.
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
            if ('abmag' in colname or 'vegamag' in colname
                    or colname in ('nn_dist', 'sharpness', 'roundness')):
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

        catalog = self.format_columns(catalog)
        catalog.meta.update(self.catalog_metadata)

        return catalog
