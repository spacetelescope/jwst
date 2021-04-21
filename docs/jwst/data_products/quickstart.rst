Quick Start Guide
=================

Transitioning from FITS (Image)
-------------------------------

If you are transitioning from using FITS to using JWST ASDF-in-FITS,
instead of doing this:

.. doctest-skip::

    >>> from astropy.io import fits
    >>> from astropy.wcs import WCS
    >>> filename = 'myfits.fits'
    >>> with fits.open(filename) as pf:
    ...     myhdu = pf['SCI']
    ...     mydata = myhdu.data
    ...     mywcs = WCS(myhdu.header)

Now, you would do this instead:

.. doctest-skip::

    >>> from jwst import datamodels
    >>> filename = 'myasdfinfits.fits'
    >>> with datamodels.open(filename) as dm:
    ...     mydata = dm.data  # Same as SCI
    ...     mywcs = dm.meta.wcs  # This is GWCS, not FITS WCS
