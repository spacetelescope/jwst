
Description
===========

jwst.assign_wcs is one of the first steps in the level 2B JWST pipeline.
It associates a WCS object with each science exposure. The WCS transforms
positions on the detector to a world coordinate frame - ICRS and wavelength.
In general there may be intermediate coordinate frames depending on the instrument.
The WCS is saved in the FITS file. It can be accessed as an attribute of
the meta object when the fits file is opened as a data model.

The forward direction of the transforms is from detector to world coordinates
and the input positions are 0-based.

The basic WCS keywords are in the primary header and the distortion
and spectral models are stored in reference files in the
`ASDF <http://asdf-standard.readthedocs.org/en/latest/>`__  format.

For each observing mode, determined by the value of ``EXP_TYPE`` in the science header,
assign_wcs retrieves reference files from CRDS and creates a pipeline of transforms from
input frame ``detector`` to a frame ``v2v3``. This part of the WCS pipeline may include
intermediate coordinate frames. The basic WCS keywords are used to create
the transform from frame ``v2v3`` to frame ``world``.

Basic WCS keywords and the transform from ``v2v3`` to ``world``
---------------------------------------------------------------

All JWST instruments use the following FITS header keywords to
define the transform from ``v2v3`` to ``world``:

``RA_REF``, ``DEC_REF`` - a fiducial point on the sky, ICRS, [deg]

``V2_REF``, ``V3_REF`` - a point in the V2V3 system which maps to ``RA_REF``, ``DEC_REF``, [arcsec]

``ROLL_REF`` - local roll angle associated with each aperture, [deg]

These quantities are used to create a 3D Euler angle rotation between the V2V3 spherical system,
associated with the telescope, and a standard celestial system.


Using the WCS interactively
---------------------------

Once a FITS file is opened as a `DataModel` the WCS can be accessed as an attribute of the meta object.
Calling it as a function with detector positions as inputs returns the
corresponding world coordinates. Using MIRI LRS fixed slit as an example:

>>> from jwst.datamodels import ImageModel
>>> exp = ImagModel(miri_fixed_assign_wcs.fits')
>>> ra, dec, lam = exp.meta.wcs(x, y)
>>> print(ra, dec, lam)
    (329.97260532549336, 372.0242999250267, 5.4176100046836675)

The WCS provides access to intermediate coordinate frames
and transforms between any two frames in the WCS pipeline in forward or
backward direction. For example, for a NIRSPEC fixed slits exposure,
which has been through the extract_2d step:

>>> exp = models.MultiSlitModel('nrs1_fixed_assign_wcs_extract_2d.fits')
>>> exp.slits[0].meta.wcs.available_frames
    ['detector', 'sca', 'bgwa', 'slit_frame', 'msa_frame', 'ote', 'v2v3', 'world']
>>> msa2detector = exp.slits[0].meta.wcs.get_transform('msa_frame', 'detector')
>>> msa2detector(0, 0, 2*10**-6)
    (5042.064255529629, 1119.8937888372516)

For each exposure, assign_wcs uses reference files and WCS header keywords
to create the WCS object. What reference files are retrieved
from CRDS is determined based on EXP_TYPE and other keywords in the science file header.


The assign_wcs step can accept any type of DataModel as input. In particular, for
multiple-integration datasets the step will accept either of these data products:
the slope results for each integration in the exposure, or the single slope image
that is the result of averaging over all integrations.

jwst.assign_wcsis based on gwcs and uses the modeling, units and coordinates subpackages in astropy.

Software dependencies:

- `gwcs <https://github.com/spacetelescope/gwcs>`__ 0.7

- `numpy <http://www.numpy.org/>`__ 1.9 or later

- `astropy <http://www.astropy.org/>`__ 1.2.1 or later

- `asdf <http://pyasdf.readthedocs.org/en/latest/>`__ 1.1.1 or later
