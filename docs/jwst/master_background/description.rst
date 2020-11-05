Description
===========
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

Inputs
------
The primary driver of the master background step is a "spec3" type Association (ASN) file
or a ``ModelContainer`` data model populated from a "spec3" ASN file. This is the same ASN file used
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
"BACKGROUND" column. The ``master_background`` step recognizes when it's working with nodded
exposures and in that case uses the data from the "BACKGROUND" column of each background
:ref:`x1d <x1d>` product to create the 1-D master background spectrum.

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
extraction region. The extracted spectral data are stored in the "FLUX" and "SURF_BRIGHT"
columns of the resulting :ref:`x1d <x1d>` product, with the "BACKGROUND" column left blank.
The ``master_background`` step recognizes when it's working with a background exposure, which is
always treated as an extended source,
and in that case uses the data from the "SURF_BRIGHT" column of each background
:ref:`x1d <x1d>` product to construct the master background spectrum.

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

NIRSpec MOS with Background Slits
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
NIRSpec MOS exposures that have one or more slits defined as background require
unique processing. The background slits take the place of both nodded and off-target
background exposures, because the background can be measured directly from the spectra
contained in those slits. Because the background spectra come from the same exposure
as the source spectra, the master background processing is applied within the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline as each individual exposure is
calibrated, rather than during :ref:`calwebb_spec3 <calwebb_spec3>` processing.
In this scenario, all source and background slits are first partially calibrated up through
the :ref:`extract_2d <extract_2d_step>` and :ref:`srctype <srctype_step>` steps of
:ref:`calwebb_spec2 <calwebb_spec2>` to produce 2D cutouts for each slit. At this
point the `master_background` step is applied, which completes the remaining calibration
steps for all of the background slits, resulting in multiple 1D extracted background
spectra. The background spectra are combined (as with other observing modes) to
create a 1D master background spectrum, and then the master background spectrum is
interpolated back into the 2D space of each source slit and subtracted. The
background-subtracted source slits then have all of their remaining
:ref:`calwebb_spec2 <calwebb_spec2>` calibration steps applied.

Note that special corrections are applied to the 2D master background data before
being subtracted from each source slit, as explained in detail below.

Creating the 1-D Master Background Spectrum
-------------------------------------------
The 1-D master background spectrum is created by combining data contained in the
:ref:`x1d <x1d>` products listed in the input ASN as ``"exptype": "background"`` members.
As noted above, the background members can be exposures of dedicated background targets
or can be a collection of exposures of a point-like source observed in a nod pattern.

For the case of dedicated background target exposures, the spectrum contained in the
"SURF_BRIGHT" column of the background :ref:`x1d <x1d>` products will be used for creating the
master background spectrum. For the case of nodded exposures, the spectrum contained
in the "BACKGROUND" column of the :ref:`x1d <x1d>` products will be used. The data
in both columns are in units of surface brightness, which is appropriate for
eventually computing the 2-D background signal.

When all the input background spectra have been collected, they are combined using the
:ref:`combine_1d <combine_1d_step>` step to produce the 1-D master background spectrum.
See the :ref:`combine_1d <combine_1d_step>` step for more details on the processes used
to create the combined spectrum.

Subtracting the Master Background
---------------------------------
The 1-D master background spectrum is interpolated by wavelength at each pixel of a 2-D source
spectrum and subtracted from it. The source data instances can be, for example, a set
of NIRSpec or MIRI IFU exposures, a set of NIRSpec MOS or fixed-slit 2-D extractions, or a set of
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

NIRSpec Background Corrections
------------------------------
Once the 1-D master background spectrum has been interpolated to the 2-D space of
the science data, NIRSpec data can sometimes need additional corrections to make the
computed background match the science data. This is due to two primary effects of NIRSpec
calibration:

- Point sources in MOS and fixed-slit mode receive wavelength offset
  corrections if the source is not centered (in the dispersion direction) within the slit.
  Hence the wavelength grid assigned to the 2-D slit cutout is shifted slightly relative
  to the wavelengths of the background signal contained in the same cutout. And because the
  flat-field, pathloss, and photom corrections/calibrations are wavelength-dependent, the
  pixel-level calibrations for the source signal are slightly different than the background.

- Point sources and uniform sources receive different pathloss and bar shadow corrections
  (in fact point sources don't receive any bar shadow correction). So the background signal
  contained within a calibrated point source cutout has received the wrong pathloss
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
For the NIRSpec IFU mode, the only effect that needs to be accounted for is the difference
between point source and uniform source pathloss corrections, because no wavelength or
bar shadow corrections are applied to IFU data. The 2-D background generated from the
master background spectrum must have the uniform source pathloss correction removed from
it and the point source pathloss correction applied to it, so the operation performed on
the IFU 2-D background is:

.. math::
 bkg(corr) = bkg * pathloss(uniform) / pathloss(point)

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

So if the primary slit (as given by the FXD_SLIT keyword) contains a point source
(as given by the SRCTYPE keyword) the corrections applied to the 2-D master background
for that slit are:

.. math::
 bkg(corr) = bkg &* [flatfield(uniform) / flatfield(point)]\\
                 &* [pathloss(uniform) / pathloss(point)]\\
                 &* [photom(point) / photom(uniform)]

For secondary slits that contain a point source, the corrections applied to the
2-D master background are simply:

.. math::
 bkg(corr) = bkg * pathloss(uniform) / pathloss(point)

NIRSpec MOS Mode
^^^^^^^^^^^^^^^^
Because the master background is subtracted from only partially calibrated slit
data, the number of correction terms is reduced relative to that for fixed slits.
At the time the 2D master background is subtracted from each source slit, those
source slits have not yet received any calibrations that are specific to point
sources and hence the background signal in the source slits has effects in it
characteristic of uniform source data. Therefore the fully-calibrated 2D
master background signal needs to have those uniform source effects imposed on it
(e.g. impose the uniform source flat-field pattern or the bar shadow pattern
into the background data),
so that it matches the actual background signal in each slit before being
subtracted. The corrections applied to the 2D master background for MOS slits are:

.. math::
 bkg(corr) = bkg &* flatfield(uniform) * pathloss(uniform)\\
                 &* barshadow(uniform) / photom(uniform)

