Reference Files
===============
The photom step uses a photom reference file and a pixel area map reference
file. The pixel area map reference file is only used when processing
imaging-mode observations.

PHOTOM CRDS Selection Criteria
------------------------------
Selection criteria for photom reference files varies a bit from instrument
to instrument:

* FGS: Instrument and Detector
* MIRI: Instrument and Detector
* NIRCam: Instrument and Detector
* NIRISS: Instrument and Detector
* NIRSpec: Instrument and EXP_TYPE

A row of data within the table that matches the mode of the science exposure
is selected by the photom step based on criteria that are instrument mode
dependent. The current row selection criteria are:

* FGS: No selection criteria (table contains a single row)
* MIRI:
   - Imager: Filter and Subarray
   - IFUs: Band
* NIRCam: Filter and Pupil
* NIRISS: Filter, Pupil, and Order number
* NIRSpec:
   - Fixed Slits: Filter, Grating, and Slit name
   - IFU and MSA: Filter and Grating

PHOTOM Reference File Format
----------------------------
Photom reference files are FITS format with a single BINTABLE extension.  The
primary data unit is always empty.  The columns of the table vary with 
instrument according to the selection criteria listed above. The first few
columns always correspond to the selection criteria, such as Filter and
Pupil, or Filter and Grating. The remaining columns contain the data relevant
to the photometric conversion and consist of PHOTMJSR, UNCERTAINTY, NELEM,
WAVELENGTH, and RELRESPONSE.

* FILTER (string) - MIRI, NIRCam, NIRISS, NIRSpec
* PUPIL (string) - NIRCam, NIRISS
* ORDER (integer) - NIRISS
* GRATING (string) - NIRSpec
* SLIT (string) - NIRSpec Fixed-Slit
* SUBARRAY (string) - MIRI Imager/LRS
* BAND (string) - MIRI MRS
* PHOTMJSR (float) - all instruments
* UNCERTAINTY (float) - all instruments
* NELEM (int) - if NELEM > 0, then NELEM entries are read from each of the
  WAVELENGTH and RELRESPONSE arrays
* WAVELENGTH (float 1-D array)
* RELRESPONSE (float 1-D array)

The primary header of the photom reference file contains the keywords PIXAR_SR
and PIXAR_A2, which give the average pixel area in units of steradians and
square arcseconds, respectively.

AREA CRDS Selection Criteria
------------------------------
Selection criteria for AREA reference files vary a bit from instrument to
instrument:

* FGS: Instrument and Detector
* MIRI: Instrument, Detector, and Exp_Type
* NIRCam: Instrument, Detector, and Exp_Type
* NIRISS: Instrument, Detector, and Exp_Type
* NIRSpec: Instrument, Detector, and Exp_Type

AREA Reference File Format
-----------------------------
Pixel area map reference files are FITS format with a single image extension
with 'EXTNAME=SCI', which contains a 2-D floating-point array of values. The FITS
primary data array is always empty. The primary header contains the keywords
PIXAR_SR and PIXAR_A2, which should have the same values as the keywords in
the header of the corresponding photom reference file.

Constructing a PHOTOM Reference File
------------------------------------
The most straight-forward way to construct a PHOTOM reference file is to
populate a photom data model within python and then save the data model to a
FITS file. Each instrument has its own photom data model, which contains the
columns of information unique to that instrument:

* NircamPhotomModel
* NirissPhotomModel
* NirspecPhotomModel
* MiriImgPhotomModel
* MiriMrsPhotomModel

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
 >>> data=np.array(zip(filter,pupil,photf,uncer,nelem,wave,resp),dtype=output.phot_table.dtype)
 >>> output.phot_table=data
 >>> output.save('niriss_photom_0001.fits')

