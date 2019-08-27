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
set to 'COMPLETE', the master background step is skipped. If the :ref:`calwebb_spec2 <calwebb_spec2>`
background step was not applied, the master background step will proceed.
The user can override this logic, if desired, by setting the step argument ``--force_subtract``
to ``True``, in which case master
background subtraction will be applied regardless of the value of S_BKDSUB (see
:ref:`msb_step_args`).

Upon successful completion of the step, the S_MSBSUB keyword is set to 'COMPLETE' in the
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

There are two main observing scenarios that are supported by this step: nodded point sources
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
source type to "POINT".

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
"BACKGROUND" column. The master_background step recognizes when it's working with nodded
exposures and in that case uses the data from the "BACKGROUND" column of each background
:ref:`x1d <x1d>` product.

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
"BKGDTARG" set to True in the header. During :ref:`calwebb_spec2 <calwebb_spec2>` processing,
the :ref:`srctype <srctype_step>` step recognizes these and sets their source type to
"EXTENDED", because all dedicated background exposures are to be processed as extended sources.

This in turn causes the :ref:`extract_1d <extract_1d_step>` step at
the end of :ref:`calwebb_spec2 <calwebb_spec2>` to extract a spectrum in extended source mode,
which uses the entire field-of-view (whether it be a slit image or an IFU cube) as the
extraction region. The extracted spectral data are stored in the "FLUX" and "SURF_BRIGHT"
columns of the resulting :ref:`x1d <x1d>` product, with the "BACKGROUND" column left blank.
The master_background step recognizes when it's working with a background exposure, which is
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
resulting from those slits. In this observing scenario, all source and background
slits are processed through all of the :ref:`calwebb_spec2 <calwebb_spec2>` steps,
resulting in extracted spectra for the sources and backgrounds. The master_background
step then combines the background spectra into a master background spectrum and
performs background subtraction on all of the source slit data, the same as the
other observing scenarios.

There are, however, additional factors that must be taken into account in order to
yield scientifically-correct results. This is due to the fact that background slits
are processed as extended sources in :ref:`calwebb_spec2 <calwebb_spec2>`, but
the background spectra must then be applied to some slits that were processed as
point sources. Extended sources and point sources receive different corrections
in steps like :ref:`pathloss <pathloss_step>` and :ref:`barshadow <barshadow_step>`,
hence the data are not a one-to-one match. This requires additional operations to
be performed on the background spectrum before it can be correctly applied to slits
containing point sources.

**These unique capabilities are not yet implemented in the master_background step,
hence NIRSpec MOS observations can not be processed properly at this time.** The
only workable option for NIRSpec MOS data at this time is to employ the
user-supplied background spectrum option, which is then subtracted from every
slit instance in a MOS exposure.

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
