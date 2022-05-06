Description
===========

:Class: `jwst.assign_mtwcs.AssignMTWcsStep`
:Alias: assign_mtwcs


The ``jwst.assign_mtwcs`` step modifies the WCS output frame in each exposure of
a Moving Target (MT) observation association, such that the WCS is centered at the
average location of the target within the whole association.
This results in proper alignment of multiple exposures, which takes place downstream
in the calibration pipeline, in the target frame, rather than the background
sky frame.

A moving target will naturally be located at different sky coordinates (RA, Dec)
across multiple exposures within an MT observation. When multiple images or
spectra get combined during Stage 3 processing, the relative alignment of the
images/spectra is based on the sky coordinates of each exposure.
In the case of moving targets, where the RA/Dec of the target is changing
between exposures, the normal alignment process would result in the target being
at different image coordinates and hence coming out either smeared (for slowly
moving targets) or at multiple locations within the combined data. This step
modifies the WCS of each exposure to recenter it at a common RA/Dec for the
target, so that subsequent image alignment and combination has the target
properly aligned.

The step is executed at the beginning of the :ref:`calwebb_image3 <calwebb_image3>`
and :ref:`calwebb_spec3 <calwebb_spec3>` pipelines, so that all subsequent steps
that rely on WCS information use the frame centered on the target.

This step depends on keywords that are unique to MT exposures, as shown in the
following table.

+------------+-------------------------+------------+-------------------------------+
| | FITS     | Data model attribute    | | Type     | Description                   |
| | Keyword  |                         | | (Value)  |                               |
+============+=========================+============+===============================+
| | TARGTYPE | | meta.target.type      | | string   | Type of target                |
|            |                         | | (moving) |                               |
+------------+-------------------------+------------+-------------------------------+
| | MT_RA    | | meta.wcsinfo.mt_ra    | | number   | | Target RA and Dec at        |
| | MT_DEC   | | meta.wcsinfo.mt_dec   | | number   | | mid-point of exposure [deg] |
+------------+-------------------------+------------+-------------------------------+
| | MT_AVRA  | | meta.wcsinfo.mt_avra  | | number   | | Target RA and Dec averaged  |
| | MT_AVDEC | | meta.wcsinfo.mt_avdec | | number   | | between exposures [deg]     |
+------------+-------------------------+------------+-------------------------------+

A "TARGTYPE" value of "moving" is used to identify exposures as containing
a moving target. The keywords "MT_RA" and "MT_DEC" are populated in the uncalibrated
(:ref:`uncal <uncal>`) product for each exposure and give the position of the
target at the mid-point of each exposure. The ``assign_mtwcs`` step computes the
average of the "MT_RA" and "MT_DEC" values across all expsoures in an association
and stores the result in the "MT_AVRA" and "MT_AVDEC" keywords of each exposure.

In addition to populating the "MT_AVRA" and "MT_AVDEC" keywords, this step adds
another transform to the original WCS in each exposure that results in the WCS
frame being centered at "MT_AVRA" and "MT_AVDEC".
The transform of the original WCS associated with the science aperture pointing
(i.e. without the additional MT correction) can be accessed by executing::

  sci_transform = model.meta.wcs.get_transform('detector', 'world')
