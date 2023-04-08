"""
`skystatistics` module provides statistics computation class used by
:py:func:`~jwst.skymatch.skymatch.match`
and :py:class:`~jwst.skymatch.skyimage.SkyImage`.

:Authors: Mihai Cara (contact: help@stsci.edu)


"""
# THIRD PARTY
from stsci.imagestats import ImageStats
from copy import deepcopy

__all__ = ['SkyStats']
__taskname__ = 'skystatistics'
__author__ = 'Mihai Cara'


class SkyStats():
    """
    This is a superclass build on top of
    :py:class:`stsci.imagestats.ImageStats`. Compared to
    :py:class:`stsci.imagestats.ImageStats`, `SkyStats` has
    "persistent settings" in the sense that object's parameters need to be
    set once and these settings will be applied to all subsequent
    computations on different data.

    """

    def __init__(self, skystat='mean', lower=None, upper=None,
                 nclip=5, lsig=4.0, usig=4.0, binwidth=0.1, **kwargs):
        """Initializes the SkyStats object.

        Parameters
        -----------
        skystat : {'mode', 'median', 'mode', 'midpt'}, optional
            Sets the statistics that will be returned by `~SkyStats.calc_sky`.
            The following statistics are supported: 'mean', 'mode', 'midpt',
            and 'median'. First three statistics have the same meaning as in
            `stsdas.toolbox.imgtools.gstatistics <http://stsdas.stsci.edu/\
cgi-bin/gethelp.cgi?gstatistics>`_
            while 'median' will compute the median of the distribution.

        lower : float, None, optional
            Lower limit of usable pixel values for computing the sky.
            This value should be specified in the units of the input image(s).

        upper : float, None, optional
            Upper limit of usable pixel values for computing the sky.
            This value should be specified in the units of the input image(s).

        nclip : int, optional
            A non-negative number of clipping iterations to use when computing
            the sky value.

        lsig : float, optional
            Lower clipping limit, in sigma, used when computing the sky value.

        usig : float, optional
            Upper clipping limit, in sigma, used when computing the sky value.

        binwidth : float, optional
            Bin width, in sigma, used to sample the distribution of pixel
            brightness values in order to compute the sky background
            statistics.

        kwargs : dict
            A dictionary of optional arguments to be passed to `ImageStats`.

        """
        self.npix = None
        self.skyval = None

        self._fields = ','.join(['npix', skystat])

        self._kwargs = deepcopy(kwargs)
        if 'fields' in self._kwargs:
            del self._kwargs['fields']
        if 'image' in self._kwargs:
            del self._kwargs['image']
        self._kwargs['lower'] = lower
        self._kwargs['upper'] = upper
        self._kwargs['nclip'] = nclip
        self._kwargs['lsig'] = lsig
        self._kwargs['usig'] = usig
        self._kwargs['binwidth'] = binwidth

        self._skystat = {'mean': self._extract_mean,
                         'mode': self._extract_mode,
                         'median': self._extract_median,
                         'midpt': self._extract_midpt
                         }[skystat]

    def _extract_mean(self, imstat):
        return imstat.mean

    def _extract_median(self, imstat):
        return imstat.median

    def _extract_mode(self, imstat):
        return imstat.mode

    def _extract_midpt(self, imstat):
        return imstat.midpt

    def calc_sky(self, data):
        """ Computes statistics on data.

        Parameters
        -----------
        data : numpy.ndarray
            A numpy array of values for which the statistics needs to be computed.

        Returns
        --------
        statistics : tuple
            A tuple of two values: (`skyvalue`, `npix`), where `skyvalue` is the statistics
            specified by the `skystat` parameter during the initialization
            of the `SkyStats` object and `npix` is the number of pixels used
            in computing the statistics reported in `skyvalue`.

        """
        imstat = ImageStats(image=data, fields=self._fields,
                            **(self._kwargs))
        self.skyval = self._skystat(imstat)
        self.npix = imstat.npix
        return self.skyval, self.npix

    def __call__(self, data):
        return self.calc_sky(data)
