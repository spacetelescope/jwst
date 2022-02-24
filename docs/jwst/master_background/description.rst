Description
===========

:Classes: `jwst.master_background.MasterBackgroundStep`, `jwst.master_background.MasterBackgroundMosStep`
:Aliases: master_background, master_background_mos

Master background subtraction is one form of background subtraction available for
spectroscopic data. See :ref:`Background Subtraction <background_subtraction>` for an
overview of all the available methods and where they occur within the various stages
of the calibration pipeline.

The master background subtraction step subtracts background signal from
2-D spectroscopic data using a 1-D master background spectrum. The 1-D master background
spectrum is created from one or more input exposures, or can alternatively be supplied
by the user. The 1-D background spectrum - surface brightness
versus wavelength - is projected into the
2-D space of source data based on the wavelength of each pixel in the 2-D data. The resulting
2-D background signal is then subtracted directly from the 2-D source data.

Logic built into the step checks to see if the exposure-based :ref:`background <background_step>`
subtraction step in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline has already been
performed on the input images, based on the value of the S_BKDSUB keyword. If S_BKGSUB is
set to "COMPLETE", the master background step is skipped. If the :ref:`calwebb_spec2 <calwebb_spec2>`
background step was not applied, the master background step will proceed.
The user can override this logic, if desired, by setting the step argument ``--force_subtract``
to ``True``, in which case master background subtraction will be applied regardless of the
value of S_BKDSUB (see :ref:`msb_step_args`).

Upon successful completion of the step, the S_MSBSUB keyword is set to "COMPLETE" in the
output product. The background-subtracted results are returned as a new data model, leaving
the input model unchanged.

Note: The application of master background subtraction to NIRSpec Fixed-Slit, IFU, and MOS
observations requires special handling, due to unique types of calibrations that are
applied to these modes. NIRSpec MOS mode requires even more special handling than NIRSpec
Fixed-Slit and IFU. The next several sections pertain primarily to MIRI MRS and LRS Fixed-Slit,
and in a general way to NIRSpec Fixed-Slit and IFU modes. Details regarding all NIRSpec
modes are given later in :ref:`NIRSpec Master Background Subtraction <nirspec_modes>`.

Inputs
------
The primary driver of the master background step is usually a `spec3` type Association (ASN) file
or a ``ModelContainer`` data model populated from a `spec3` ASN file. This is the same ASN file used
as input to the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline, which defines a stage 3 combined product
and its input members. The list of input members includes both "science" and "background"
exposure types. The master background subtraction step uses the input members designated
with ``"exptype": "background"`` to create the master background spectrum (see example_asn1_).
These need to be :ref:`x1d <x1d>` products created from individual exposures at the end of
the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline, containing spectra of background regions.
The master background signal will be subtracted from all input members designated as
``"exptype": "science"`` in the ASN, resulting in a new version of each science input. These inputs
need to be :ref:`cal <cal>` products created from individual exposures by the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

There are two main observing scenarios that are supported by this step: nodded exposures of point sources
and off-source background exposures of extended targets. A third type of operation is performed
for NIRSpec MOS observations that include background slits. The details for each mode are explained
below.

Nodded Point Sources
^^^^^^^^^^^^^^^^^^^^
If an observation uses a nodding type dither pattern to move a small or point-like source within
the field-of-view, it is assumed that part of the field-of-view in each exposure is also suitable
for measuring background. Exposures of this type are identified by the pipeline based on their
"PATTTYPE" (primary dither pattern type) keyword value. The value will either contain the
substring "NOD" somewhere within the name (e.g. "2-POINT-NOD" or "ALONG-SLIT-NOD"), or will
be set to "POINT-SOURCE" (for MIRI MRS).  The :ref:`calwebb_spec2 <calwebb_spec2>`
:ref:`srctype <srctype_step>` step recognizes these PATTTYPE values and sets the
source type to "POINT."

This in turn causes the :ref:`extract_1d <extract_1d_step>` step at
the end of :ref:`calwebb_spec2 <calwebb_spec2>` to extract spectra for both source and
background regions. For IFU exposures the background region is typically an annulus that is
concentric with a circular source region. For slit-like modes, one or more background regions can
be defined in the :ref:`extract1d <extract1d_reffile>` reference file, flanking the central source region.
In both cases, the extraction regions are centered within
the image/cube at the RA/Dec of the target. Hence for nodded exposures, the location of the
extraction regions follows the movement of the source in each exposure. The extracted
data from the source region are stored in the "FLUX" and "SURF_BRIGHT" (surface brightness)
columns of the :ref:`x1d <x1d>` product, while the background extraction is stored in the
"BACKGROUND" column. The ``master_background`` step uses the data from the "BACKGROUND" column
of each background :ref:`x1d <x1d>` product to create the 1-D master background spectrum.

Below is an example ASN file for a simple 2-point nodded observation consisting of two
exposures.

.. _example_asn1:

::

  {
      "asn_type": "spec3",
      "asn_rule": "candidate_Asn_IFU",
      "program": "00626",
      "asn_id": "c1003",
      "target": "t001",
      "asn_pool": "jw00626_20190128T194403_pool",
      "products": [
          {"name": "jw00626-c1003_t001_nrs",
              "members": [
                  {"expname": "jw00626009001_02101_00001_nrs1_cal.fits",
                    "exptype": "science",
                    "asn_candidate": "('c1003', 'background')"
                  },
                  {"expname": "jw00626009001_02102_00001_nrs1_cal.fits",
                   "exptype": "science", 
                   "asn_candidate": "('c1003', 'background')"
                  },
                  {"expname": "jw00626009001_02101_00001_nrs1_x1d.fits",
                   "exptype": "background",
                   "asn_candidate": "('c1003', 'background')"
                  },
                  {"expname": "jw00626009001_02102_00001_nrs1_x1d.fits",
                   "exptype": "background",
                   "asn_candidate": "('c1003', 'background')"
                  }
              ]
          }
      ]
  }

As you can see, the same two exposures are defined as
being both "science" and "background" members, because they both contain the target
of interest and a region of background. The "science" members, which are the
:ref:`cal <cal>` products created by the :ref:`calwebb_spec2 <calwebb_spec2>`
pipeline, are the data files that will have the master background subtraction
applied, while the "background" members are the :ref:`x1d <x1d>` spectral
products from which the master background spectrum will be created.
The combined master background spectrum will be subtracted from each of the 
two science exposures.

Extended Source with Dedicated Background Exposures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Observations of extended sources must obtain exposures of a separate background target/field in
order to measure the background. Exposures of a background target are identified by the keyword
"BKGDTARG" set to `True` in the header. During :ref:`calwebb_spec2 <calwebb_spec2>` processing,
the :ref:`srctype <srctype_step>` step recognizes these and sets their source type to
"EXTENDED", because all dedicated background exposures are to be processed as extended sources.

This in turn causes the :ref:`extract_1d <extract_1d_step>` step at
the end of :ref:`calwebb_spec2 <calwebb_spec2>` to extract a spectrum in extended source mode,
which uses the entire field-of-view (whether it be a slit image or an IFU cube) as the
extraction region.
The ``master_background`` step recognizes which type of background exposure it's working with
and uses the appropriate data from the :ref:`x1d <x1d>` product to construct the master
background spectrum.

Below is an example ASN file for an extended source observation that includes background target
exposures, using a 2-point dither for both the science and background targets.

.. _example_asn2:

::

  {
      "asn_type": "spec3",
      "asn_rule": "candidate_Asn_IFU",
      "program": "00626",
      "asn_id": "c1004",
      "target": "t002",
      "asn_pool": "jw00626_20190128T194403_pool",
      "products": [
          {"name": "jw00626-c1004_t002_nrs",
              "members": [
                  {"expname": "jw00626009001_02101_00001_nrs1_cal.fits",
                    "exptype": "science",
                    "asn_candidate": "('c1004', 'background')"
                  },
                  {"expname": "jw00626009001_02102_00001_nrs1_cal.fits",
                   "exptype": "science", 
                   "asn_candidate": "('c1004', 'background')"
                  },
                  {"expname": "jw00626009001_02103_00001_nrs1_x1d.fits",
                   "exptype": "background",
                   "asn_candidate": "('c1004', 'background')"
                  },
                  {"expname": "jw00626009001_02104_00001_nrs1_x1d.fits",
                   "exptype": "background",
                   "asn_candidate": "('c1004', 'background')"
                  }
              ]
          }
      ]
  }

In this example there are two exposures of the science target, labeled as "science"
members, and two exposures of the background target, labeled as "background"
members. As before, the science members use :ref:`cal <cal>` products as input
and the background members use :ref:`x1d <x1d>` products as input.
The master background step will first combine the data from the two background
members into a master background spectrum and then subtract it from each of the
two science exposures.

Creating the 1-D Master Background Spectrum
-------------------------------------------
The 1-D master background spectrum is created by combining data contained in the
:ref:`x1d <x1d>` products listed in the input ASN as ``"exptype": "background"`` members.
As noted above, the background members can be exposures of dedicated background targets
or can be a collection of exposures of a point-like source observed in a nod pattern.

When all of the input background spectra have been collected, they are combined using the
:ref:`combine_1d <combine_1d_step>` step to produce the 1-D master background spectrum.
See the :ref:`combine_1d <combine_1d_step>` step for more details on the processes used
to create the combined spectrum.

Subtracting the Master Background
---------------------------------
The 1-D master background spectrum is interpolated by wavelength at each pixel of a 2-D source
spectrum and subtracted from it. The source data instances can be, for example, a set
of NIRSpec or MIRI IFU exposures, a set of NIRSpec fixed-slit 2-D extractions, or a set of
nodded MIRI LRS fixed-slit exposures. The subtraction is performed on all data instances
within all input science exposures. For example, if there are 3 NIRSpec fixed-slit exposures,
each containing data from multiple slits, the subtraction is applied one-by-one to all slit
instances in all exposures. For each data instance to be subtracted the following steps are
performed:

- Compute a 2-D wavelength grid corresponding to the 2-D source data. For some observing modes,
  such as NIRSpec MOS and fixed-slit, a 2-D wavelength array is already computed and attached to the data
  in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline :ref:`extract_2d <extract_2d_step>` step.
  If such a wavelength array is present, it is used. For modes that don't have a 2-D
  wavelength array contained in the data, it is computed on the fly using the WCS object
  for each source data instance.

- Compute the background signal at each pixel in the 2-D wavelength grid by interpolating within
  the 1-D master background spectrum as a function of wavelength.
  Pixels in the 2-D source data with an undefined wavelength (e.g. wavelength array value
  of NaN) or a wavelength that is beyond the limits of the master background spectrum receive
  special handling. The interpolated background value is set to zero and a DQ flag of
  "DO_NOT_USE" is set.

- Subtract the resulting 2-D background image from the 2-D source data. DQ values from the
  2-D background image are propagated into the DQ array of the subtracted science data.

.. _nirspec_modes:

NIRSpec Master Background Corrections
-------------------------------------
The master background subtraction methods and processing flow for NIRSpec Fixed-Slit
and IFU modes is largely the same as what's outlined above, with some additional
operations that need to be applied to accommodate some of the unique calibrations
applied to NIRSpec data. NIRSpec MOS mode requires even more special handling.
This is due to two primary effects of NIRSpec calibration:

- Point sources in MOS and Fixed-Slit mode receive wavelength offset
  corrections if the source is not centered (along the dispersion direction) within the slit.
  Hence the wavelength grid assigned to each 2-D slit cutout can be shifted slightly relative
  to the wavelengths of the background signal contained in the same cutout. And because the
  flat-field, pathloss, and photom corrections/calibrations are wavelength-dependent, the
  pixel-level calibrations for the source signal are slightly different than the background.

- Point sources and uniform sources receive different pathloss and bar shadow corrections
  (in fact point sources don't receive any bar shadow correction). So the background signal
  contained within a calibrated point source cutout has received a different pathloss
  correction and hasn't received any bar shadow correction. Meanwhile, the master background
  is created from data that had corrections for a uniform source applied to it and hence
  there's a mismatch relative to the point source data.

The 2-D background that's initially created from the 1-D master background is essentially
a perfectly calibrated background signal. However, due to the effects mentioned above, the
actual background signal contained within a calibrated point source slit (or IFU image) is not
perfect (e.g. it still has the bar shadow effects in it). So all of these effects need to be
accounted for in the computed 2-D background before subtracting from the source data.

NIRSpec IFU Mode
^^^^^^^^^^^^^^^^
For the NIRSpec IFU mode, the overall processing flow is the same as other modes, in that
the 1-D master background spectrum is created and applied during
:ref:`calwebb_spec3 <calwebb_spec3>` processing, as outlined above.
No wavelength offset or bar shadow corrections are applied to IFU data, so any differences
due to the way those calibrations are applied are not relevant to IFU mode. So the only
effect that needs to be accounted for in the 2-D background generated
from the master background is the difference between point source and uniform source
pathloss corrections. This is accomplished by removing the uniform source pathloss correction
from the 2-D background signal and applying the point source pathloss correction to it. It
is then in a state where it matches the background signal contained in the point source IFU
image from which it will be subtracted.
Mathematically, the operation performed on the IFU 2-D background is:

.. math::
 bkg(corr) = bkg * pathloss(uniform) / pathloss(point)

The uniform and point source pathloss correction arrays referenced above are
retrieved from the :ref:`cal <cal>` products used as input to the master background
step. They are computed by the :ref:`pathloss <pathloss_step>` step during
:ref:`calwebb_spec2 <calwebb_spec2>` processing and stored as extra extensions in
the :ref:`cal <cal>` products.

NIRSpec Fixed-Slit Mode
^^^^^^^^^^^^^^^^^^^^^^^
NIRSpec fixed slit data receive flat-field, pathloss, and photometric calibrations,
all of which are wavelength-dependent, and the pathloss correction is also source
type dependent. Fixed slit data do not receive a bar shadow correction. Only slits
containing a point source can have a wavelength correction applied, to account for
source centering within the slit, hence slits containing uniform sources receive
the same flat-field and photometric calibrations as background spectra and
therefore don't require corrections for those two calibrations. Furthermore, the
source position in the slit is only known for the primary slit in an exposure, so
even if the secondary slits contain point sources, no wavelength correction can
be applied, and therefore again the flat-field and photometric calibrations are
the same as for background spectra. This means only the pathloss correction
difference between uniform and point sources needs to be accounted for in the
secondary slits.

Therefore if the primary slit (as given by the FXD_SLIT keyword) contains a point source
(as given by the SRCTYPE keyword) the corrections that need to be applied to the 2-D
master background for that slit are:

.. math::
 bkg(corr) = bkg &* [flatfield(uniform) / flatfield(point)]\\
                 &* [pathloss(uniform) / pathloss(point)]\\
                 &* [photom(point) / photom(uniform)]

For secondary slits that contain a point source, the correction applied to the
2-D master background is simply:

.. math::
 bkg(corr) = bkg * pathloss(uniform) / pathloss(point)

The uniform and point source versions of the flat-field, pathloss, and photom
corrections are retrieved from the input :ref:`cal <cal>` product. They
are computed and stored there during the execution of each of those steps
during :ref:`calwebb_spec2 <calwebb_spec2>` processing of NIRSpec Fixed-Slit
exposures.

NIRSpec MOS Mode
^^^^^^^^^^^^^^^^
Master background subtraction for NIRSpec MOS mode shares the high-level concepts
of other modes, but differs greatly in the details. Most importantly, the source
of the master background spectrum does not come from either nodded exposures or
exposures of a background target. The background data instead come from designated
background MSA slitlets contained with the same exposure as the science targets.
Alternatively, a user can supply a master background spectrum to be used, as is
the case for all other modes.
The master background processing for MOS mode is therefore done within the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline when processing individual MOS
exposures, rather than in the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline.
Applying the master background subtraction within the :ref:`calwebb_spec2 <calwebb_spec2>`
pipeline also has advantages due to the complex series of operations that need
to be performed, as described below.

During :ref:`calwebb_spec2 <calwebb_spec2>` processing, all source and background
slits are first partially calibrated up through the :ref:`extract_2d <extract_2d_step>`
and :ref:`srctype <srctype_step>` steps of :ref:`calwebb_spec2 <calwebb_spec2>`,
which results in 2D cutouts for each slit with the source type identified. At this
point the `master_background_mos` step is applied, which is a unique version
of the step specifically tailored to NIRSpec MOS mode. 

This version of the master background step completes the remaining calibration
for all slits, but treats them all as extended sources and saves the correction
arrays from each step (e.g. flat-field, pathloss, photom) for each slit, so that
they can be used later to apply corrections to the background data. The resulting
extracted 1D spectra from the background slits are combined to create the
master background spectrum. The master background spectrum is then interpolated
into the 2D space of each slit and has the photom, barshadow, pathloss, and
flat-field corrections removed from the 2D background arrays, so that the
background data now match the partially calibrated slit data from which they'll
be subtracted. Mathematically, the corrections applied to the 2D master background
for each MOS slit are:

.. math::
 bkg(corr) = bkg &* flatfield(uniform) * pathloss(uniform)\\
                 &* barshadow(uniform) / photom(uniform)

Once the corrected 2D backgrounds have been subtracted from each slit,
processing returns to the :ref:`calwebb_spec2 <calwebb_spec2>` flow, where all
of the remaining calibration steps are applied to each slit, resulting in
background-subtracted and fully calibrated 2D cutouts (:ref:`cal <cal>` and
:ref:`s2d <s2d>` products) and extracted 1D spectra (:ref:`x1d <x1d>` products).

The detailed list of operations performed when applying master background
subtraction to MOS data during :ref:`calwebb_spec2 <calwebb_spec2>` processing is
as follows:

1) Process all slitlets in the MOS exposure up through the
   :ref:`extract_2d <extract_2d_step>` and :ref:`srctype <srctype_step>` steps
2) The `master_background_mos` step temporarily applies remaining calibration
   steps up through :ref:`photom <photom_step>` to all slits, treating them all as
   extended sources (appropriate for background signal), and saving the extended
   source correction arrays for each slit in an internal copy of the data model
3) If a user-supplied master background spectrum is **not** given, the
   :ref:`resample_spec <resample_step>` and :ref:`extract_1d <extract_1d_step>`
   steps are applied to the calibrated background slits, resulting
   in extracted 1D background spectra
4) The 1D background spectra are combined, using the
   :ref:`combine_1d <combine_1d_step>` step, into a master background spectrum
5) If a user-supplied master background **is** given, steps 3 and 4 are skipped and
   the user-supplied spectrum is inserted into the processing flow
6) The master background spectrum (either user-supplied or created on-the-fly) is
   expanded into the 2D space of each slit
7) The 2D background "image" for each slit is processed in **inverse** mode through
   the :ref:`photom <photom_step>`, :ref:`barshadow <barshadow_step>`,
   :ref:`pathloss <pathloss_step>`, and :ref:`flatfield <flatfield_step>` steps,
   using the correction arrays that were computed in step 2, so that the background
   data now matches the partially calibrated background signal in each slit
8) The corrected 2D background is subtracted from each slit
9) The background-subtracted slits are processed through all remaining
   :ref:`calwebb_spec2 <calwebb_spec2>` calibration steps, using the corrections
   appropriate for the source type in each slit
