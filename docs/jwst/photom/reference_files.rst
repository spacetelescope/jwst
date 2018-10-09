Reference Files
===============
The photom step uses a photom reference file and a pixel area map reference
file. The pixel area map reference file is only used when processing
imaging and NIRSpec IFU observations.

.. include:: ../includes/standard_keywords.rst

.. include:: photom_selection.rst
             

PHOTOM Row Selection
^^^^^^^^^^^^^^^^^^^^

A row of data within the table that matches the mode of the science exposure
is selected by the photom step based on criteria that are instrument mode
dependent. The current row selection criteria are:

* FGS: No selection criteria (table contains a single row)
* MIRI:
   - Imager (includes LRS): Filter and Subarray
   - MRS: Does not use table-based reference file
* NIRCam: Filter and Pupil
* NIRISS:
   - Imaging: Filter and Pupil
   - Spectroscopic: Filter, Pupil, and Order number
* NIRSpec:
   - Fixed Slits: Filter, Grating, and Slit name
   - IFU and MOS: Filter and Grating

.. include:: area_selection.rst
     
Reference File Format
---------------------

PHOTOM Reference File Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Except for MIRI MRS, photom reference files are FITS format with a single
BINTABLE extension named PHOTOM.  The primary data array is always empty.  The
columns of the table vary with instrument according to the selection criteria
listed above. The first few columns always correspond to the selection
criteria, such as Filter and Pupil, or Filter and Grating. The remaining
columns contain the data relevant to the photometric conversion and consist of
PHOTMJSR, UNCERTAINTY, NELEM, WAVELENGTH, and RELRESPONSE.  The table column
names and data types are listed below.


* FILTER (string) - MIRI, NIRCam, NIRISS, NIRSpec
* PUPIL (string) - NIRCam, NIRISS
* ORDER (integer) - NIRISS
* GRATING (string) - NIRSpec
* SLIT (string) - NIRSpec Fixed-Slit
* SUBARRAY (string) - MIRI Imager/LRS
* PHOTMJSR (float) - all instruments
* UNCERTAINTY (float) - all instruments
* NELEM (int) - if NELEM > 0, then NELEM entries are read from each of the
  WAVELENGTH and RELRESPONSE arrays
* WAVELENGTH (float 1-D array)
* RELRESPONSE (float 1-D array)

The primary header of the photom reference file contains the keywords PIXAR_SR
and PIXAR_A2, which give the average pixel area in units of steradians and
square arcseconds, respectively.

MIRI MRS Photom Reference File Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MIRI MRS photom reference files do not contain tabular information,
but instead contain the following FITS extensions:

* SCI  IMAGE  2D float
* ERR  IMAGE  2D float
* DQ   IMAGE  2D unsigned-integer
* DQ_DEF  TABLE
* PIXSIZ  IMAGE  2D float

The SCI extension contains a 2D array of spectral sensitivity factors
corresponding to each pixel in a 2D MRS slice image. The sensitivity factors
are in units of DN/sec/mJy/pixel. The ERR extension contains a 2D array of
uncertainties for the SCI values, in the same units. The DQ extension
contains a 2D array of bit-encoded data quality flags for the SCI values.
The DQ_DEF extension contains a table listing the definitions of the values
used in the DQ array. The PIXSIZ extension contains a 2D array of pixel
sizes, in units of square-arcsec.

The SCI and PIXSIZ array values are both divided into the science product
SCI and ERR arrays, yielding image pixels that are units of mJy/sq-arcsec.

Scalar PHOTMJSR and PHOTUJA2 values are stored in primary header keywords
in the MIRI MRS photom reference files and are copied into the science
product header by the photom step.


AREA Reference File Format
^^^^^^^^^^^^^^^^^^^^^^^^^^
Pixel area map reference files are FITS format with a single image extension
with 'EXTNAME=SCI', which contains a 2-D floating-point array of values. The
FITS primary data array is always empty. The primary header contains the
keywords PIXAR_SR and PIXAR_A2, which should have the same values as the
keywords in the header of the corresponding photom reference file.

Constructing a PHOTOM Reference File
------------------------------------
The most straight-forward way to construct a PHOTOM reference file is to
populate a photom data model within python and then save the data model to a
FITS file. Each instrument has its own photom data model, which contains the
columns of information unique to that instrument:

* FgsPhotomModel
* NircamPhotomModel
* NirissPhotomModel
* NirspecPhotomModel (NIRSpec imaging, IFU, MOS)
* NirspecFSPhotomModel (NIRSpec fixed slits)
* MiriImgPhotomModel (MIRI imaging)
* MiriMrsPhotomModel (MIRI MRS)

A NIRISS photom reference file, for example, could be constructed as follows
from within the python environment::

 >>> from jwst import models
 >>> import numpy as np
 >>> output=models.NirissPhotomModel()
 >>> filter=np.array(['F277W','F356W','CLEAR'])
 >>> pupil=np.array(['CLEARP','CLEARP','F090W'])
 >>> photf=np.array([1.e-15,2.e-15,3.e-15])
 >>> uncer=np.array([1.e-17,2.e-17,3.e-17])
 >>> nelem=np.zeros(3)
 >>> wave=np.zeros(3)
 >>> resp=np.zeros(3)
 >>> data=np.array(list(zip(filter,pupil,photf,uncer,nelem,wave,resp)),dtype=output.phot_table.dtype)
 >>> output.phot_table=data
 >>> output.save('niriss_photom_0001.fits')

