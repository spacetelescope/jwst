Description
===========

jwst.assign_wcs is one of the first steps in the level 2B JWST pipeline.
It assigns a WCS object and associates to each exposure. 
The forward direction is from detector to world coordinates.
The WCS is saved as an attribute of the meta object of a model.
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
The instrument/mode specific sections list all keywords which
are used to determine what the rules to retrieve a reference file are for each supported mode.
The transforms defined in the reference files are chained and the WCS object is saved in the
science FITS file as part of the ASDF extension. All reference files in build 5 (except one)
are in ASDF format.

The assign_wcs step can accept any type of DataModel as input. In particular, for
multiple-integration datasets the step will accept either of these data products:
the slope results for each integration in the exposure, or the single slope image
that is the result of averaging over all integrations.

jwst.assign_wcs uses the modeling, units and coordinates subpackages in astropy.

Software dependencies:

- `gwcs <https://github.com/spacetelescope/gwcs>`__ 0.7 

- `numpy <http://www.numpy.org/>`__ 1.9 or later

- `astropy <http://www.astropy.org/>`__ 1.2.1 or later

- `asdf <http://pyasdf.readthedocs.org/en/latest/>`__ 1.1.1 or later

