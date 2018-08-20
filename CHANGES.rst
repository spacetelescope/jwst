0.11.0(Unreleased)
==================

ami
---

assign_wcs
----------

- Fixed a bug in ``get_msa_open_slits`` which prevented the code
  from finding the msa metafile.                                 [#2322]

- Fixed a bug in computing the slit_y locations for Nirspec MSA
  slitlets with more than one shutter.                           [#2325]

associations
------------

- Implemented Rule for Level 2 Nirspec Fixed Slit background. [#2307]

background
----------

barshadow
---------

- Fixed a bug in ``bar_shadow.py`` interpolate() that caused
  array index to be nan                                        [#2384]

combine_1d
----------

coron
-----

csv_tools
---------

cube_build
----------

cube_skymatch
-------------

dark_current
------------

datamodels
----------

- The ``DataModel`` ``__hasattr__`` mehtod has been replaced by ``hasattr``.
  The latter creates the attribute when it is accessed.                       [#2275]

- Improved error messaging when loading fits files into data models. [#2298]

-

dq_init
-------

emission
--------

engdblog
--------

exp_to_source
-------------

extract_1d
----------

extract_2d
----------

firstframe
----------

fits_generator
--------------

- updated pyparsing to v 2.2.0 [#2382]

flatfield
---------

fringe
------

gain_scale
----------

group_scale
-----------

guider_cds
----------

imprint
-------

ipc
---

jump
----

jwpsf
-----

lastframe
---------

lib
---

linearity
---------

model_blender
-------------

mrs_imatch
----------

msaflagopen
-----------

outlier_detection
-----------------

pathloss
--------

persistence
-----------

photom
------

pipeline
--------

- Fixed a typo in calspec2 which prevented the srctype
  step from running. [#2318]

ramp_fitting
------------

refpix
------

resample
--------

reset
-----

rscd
----

saturation
----------

skymatch
--------

source_catalog
--------------

srctype
-------

stpipe
------

straylight
----------

superbias
---------

timeconversion
--------------

transforms
----------

tso_photometry
--------------

tweakreg
--------

wfs_combine
-----------

white_light
-----------

wiimatch
--------

