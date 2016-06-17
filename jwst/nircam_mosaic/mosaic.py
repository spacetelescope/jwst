"""This module contains tools for NIRCAM image mosaic."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern.six.moves import map

# STDLIB
import os

# THIRD-PARTY
import numpy as np
from astropy.io import fits
from scipy.ndimage.interpolation import zoom

__all__ = ['NircamMosaic']


class NircamMosaic(object):
    """Class to handle NIRCAM mosaic.

    Based on IDL script and its accompanying document at
    https://confluence.stsci.edu/display/JWST/FPA+mosaic

    .. note:: Currently does not support WCS.

    Parameters
    ----------
    data_ext : int, tuple, or str
        Extension, as acceptable by :func:`astropy.io.fits.getdata`,
        of data to be mosaicked.

    sw_sca_size : int
        Output size of a single SCA in the SHORT channel
        in the mosaic. SCA of LONG channel will be
        automatically resized to twice this value.

    Examples
    --------
    >>> images = ['myimage1.fits', 'myimage2.fits', ...]
    >>> my_mosaic = NircamMosaic()
    >>> mosaiclist = my_mosaic.make_mosaic(images)

    """
    _sca_size = 2048  # Actual dimension of a detector in pixels
    _sca_gap = 187.5  # Actual gap size between SW detectors in pixels
    _module_gap = 1250.0  # Actual gap size between modules in pixels
    _pri_ext = 'PRIMARY'

    def __init__(self, data_ext=('SCI', 1), sw_sca_size=100):
        self.data_ext = data_ext
        self.sw_sca_size = sw_sca_size  # Sets multiple attributes at once

    @property
    def sw_sca_size(self):
        """Output size of a single SCA in the SHORT channel
        in the mosaic. SCA of LONG channel will be
        automatically resized to about twice this value.

        """
        return self._sw_sca_size

    @sw_sca_size.setter
    def sw_sca_size(self, val):
        """Set output SHORT SCA size and other mosaic dimensions."""
        self._sw_zoom_factor = val / self._sca_size
        self._sw_sca_size = int(val)
        self._mos_sca_gap = int(self._sca_gap * self._sw_zoom_factor)
        self._mos_module_gap = int(self._module_gap * self._sw_zoom_factor)
        self._lw_zoom_factor = ((2 * self._sw_sca_size + self._mos_sca_gap) /
                                self._sca_size)

    @property
    def sw_zoom_factor(self):
        """Zoom factor for SHORT SCA."""
        return self._sw_zoom_factor

    @property
    def lw_zoom_factor(self):
        """Zoom factor for LONG SCA."""
        return self._lw_zoom_factor

    @property
    def sca_gap(self):
        """Output size of gap between SCA detectors in the mosaic."""
        return self._mos_sca_gap

    @property
    def module_gap(self):
        """Output size of gap between modules and channels in the mosaic."""
        return self._mos_module_gap

    def _get_position(self, key):
        """Mosaic position by detector or channel."""
        key = key.upper()
        if key in ('NRCA4', 'NRCB1'):
            pos = 'lower_left'
        elif key in ('NRCA2', 'NRCB3'):
            pos = 'lower_right'
        elif key in ('NRCA3', 'NRCB2'):
            pos = 'upper_left'
        elif key in ('NRCA1', 'NRCB4'):
            pos = 'upper_right'
        elif key in ('NRCALONG', 'NRCBLONG'):
            pos = 'top'
        elif key == 'SHORT':
            pos = 'bottom'
        else:
            raise ValueError('Undefined mosaic position for {0}'.format(key))
        return pos

    def _single_sw_module(self, images):
        """Mosaic a single SHORT module."""
        mosaic = None

        for im in images:
            with fits.open(im) as pf:
                dat = pf[self.data_ext].data
                detector = pf[self._pri_ext].header.get('DETECTOR', '')

            dat = zoom(dat, self.sw_zoom_factor)

            if mosaic is None:
                mosaic = np.zeros((dat.shape[0] * 2 + self.sca_gap,
                                   dat.shape[1] * 2 + self.sca_gap))

            _insert_image(self._get_position(detector), dat, mosaic)

        return mosaic

    def _combine_sw_lw(self, sw_mosaic, lw_image):
        """Combine SHORT and LONG for a single module."""
        if lw_image is None:
            return sw_mosaic

        with fits.open(lw_image) as pf:
            lw_detector = pf[self._pri_ext].header.get('DETECTOR', '')
            lw_data = pf[self.data_ext].data

        lw_data = zoom(lw_data, self.lw_zoom_factor)

        if sw_mosaic is None:
            return lw_data

        nx = max(sw_mosaic.shape[1], lw_data.shape[1])
        ny = sw_mosaic.shape[0] + lw_data.shape[0] + self.module_gap
        mosaic = np.zeros((ny, nx))
        _insert_image(self._get_position('SHORT'), sw_mosaic, mosaic)
        _insert_image(self._get_position(lw_detector), lw_data, mosaic)

        return mosaic

    def _combine_modules(self, module_a, module_b):
        """Combine mosaics of Modules A and B."""
        if module_a is None:
            return module_b
        if module_b is None:
            return module_a

        # Fast
        if module_a.shape[0] == module_b.shape[0]:
            gap = np.zeros((module_a.shape[0], self.module_gap))
            mosaic = np.hstack((module_a, gap, module_b))

        # Slow
        else:
            nx = module_a.shape[1] + module_b.shape[1] + self.module_gap
            ny = max(module_a.shape[0], module_b.shape[0])
            mosaic = np.zeros((ny, nx))

            # Insert Module A
            mosaic[:module_a.shape[0], :module_a.shape[1]] = module_a

            # Insert Module B
            x2 = mosaic.shape[1]
            x1 = x2 - module_b.shape[1]
            mosaic[:module_b.shape[0], x1:x2] = module_b

        return mosaic

    def get_single_mosaic_array(self, images):
        """Construct mosaic from images that belong to the same dataset.

        Parameters
        ----------
        images : list
            List of filenames from the same dataset.

        Returns
        -------
        mosaic : ndarray
            Mosaic image.

        """
        # Separate Modules A and B, Channels SHORT and LONG
        mod_list = {}
        for im in images:
            hdr = fits.getheader(im, self._pri_ext)
            detector = hdr.get('DETECTOR', '').upper()  # NRC[A/B][1-4/LONG]
            try:
                module = detector[3]
                channel = detector[4:]
            except IndexError:
                continue
            if channel in ('1', '2', '3', '4'):
                channel = 'SHORT'
            if ((hdr.get('INSTRUME', '').upper() != 'NIRCAM') or
                    (channel not in ('SHORT', 'LONG')) or
                    (module not in ('A', 'B'))):
                continue
            key = (channel, module)
            if key not in mod_list:
                mod_list[key] = [im]
            else:
                mod_list[key].append(im)

        # Module A mosaic
        mos_a = self._combine_sw_lw(
            self._single_sw_module(mod_list.get(('SHORT', 'A'), [])),
            mod_list.get(('LONG', 'A'), [None])[0])

        # Module B mosaic
        mos_b = self._combine_sw_lw(
            self._single_sw_module(mod_list.get(('SHORT', 'B'), [])),
            mod_list.get(('LONG', 'B'), [None])[0])

        return self._combine_modules(mos_a, mos_b)

    def make_mosaic(self, images, outpath='', outsuffix='mosaic',
                    clobber=False, debug=False):
        """Construct one mosaic for each dataset, for multiple datasets.

        Images are sorted into datasets by JWST naming convention,
        ``jw<PPPPP><OOO><VVV>_<GGSAA>_<EEEEE>_<detector>_<suffix>.fits``,
        where the ROOTNAME is defined as
        ``jw<PPPPP><OOO><VVV>_<GGSAA>_<EEEEE>``.
        Each mosaic is saved as ``ROOTNAME_<outsuffix>.fits``,
        a single-extension FITS image.

        Parameters
        ----------
        images : list
            List of filenames.

        outpath : str
            Output directory. If not given, it is the current
            working directory.

        outsuffix : str
            Output suffix.

        clobber : bool
            If `True`, overwrite existing mosaic file(s).

        debug : bool
            If `True`, print extra information to screen.

        Returns
        -------
        mosaiclist : list
            List of mosaic filenames.

        """
        # Separate different datasets
        root_list = {}
        for im in images:
            rootname = os.path.basename('_'.join(im.split('_')[:-2]))
            if rootname not in root_list:
                root_list[rootname] = [im]
            else:
                root_list[rootname].append(im)

        # Process each dataset
        def _mosaic_one(rootname):
            imlist = root_list[rootname]
            outname = os.path.join(
                outpath, '{0}_{1}.fits'.format(rootname, outsuffix))

            # Avoid regenerating mosaic if already exist.
            # This also avoids crashing at the very end.
            if not clobber and os.path.exists(outname):
                if debug:
                    print('Using existing {0}'.format(outname))
                return outname

            mosaic = self.get_single_mosaic_array(imlist)
            if mosaic is None:
                if debug:
                    print('No mosaic for {0}'.format(imlist))
                return ''

            hdu = fits.PrimaryHDU(mosaic)

            # Inherit some keywords from primary header from first image in list
            prihdr = fits.getheader(imlist[0])
            for key in ('ROOTNAME', 'TARGNAME', 'INSTRUME',
                        'FILTER', 'PUPIL', 'DATE-OBS', 'TIME-OBS'):
                if key not in prihdr:
                    continue
                hdu.header[key] = prihdr[key]

            hdu.header.add_history('Mosaic from {0}'.format(','.join(imlist)))
            hdu.writeto(outname, clobber=clobber)
            return outname

        mosaiclist = sorted(map(_mosaic_one, list(root_list.keys())))
        return [m for m in mosaiclist if m]


def _insert_image(position, dat, mosaic):
    """Insert data into appropriate position in mosaic,
    which is modified in-place.

    """
    if position in ('lower_left', 'bottom', 'left'):
        x1 = 0
        x2 = dat.shape[1]
        y1 = 0
        y2 = dat.shape[0]
    elif position in ('lower_right', 'right'):
        x2 = mosaic.shape[1]
        x1 = x2 - dat.shape[1]
        y1 = 0
        y2 = dat.shape[0]
    elif position in ('upper_left', 'top'):
        x1 = 0
        x2 = dat.shape[1]
        y2 = mosaic.shape[0]
        y1 = y2 - dat.shape[0]
    elif position in ('upper_right',):
        x2 = mosaic.shape[1]
        x1 = x2 - dat.shape[1]
        y2 = mosaic.shape[0]
        y1 = y2 - dat.shape[0]
    else:
        raise ValueError('Invalid position ({0})'.format(position))

    mosaic[y1:y2, x1:x2] = dat
