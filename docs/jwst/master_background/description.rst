Description
===========
The master background subtraction step subtracts background signal from
2-D spectroscopic data using a 1-D master background spectrum. The 1-D master background
spectrum is computed from one or more input exposures, or can alternatively be supplied
by the user. The 1-D background spectrum - flux versus wavelength - is projected into the
2-D space of source data based on the wavelength of each pixel in the 2-D data. The resulting
2-D background signal is then subtracted directly from the 2-D source data.

Upon successful completion of the step, the S_MSBSUB keyword is set to 'COMPLETE' in the
output product. The background-subtracted results are returned as a new data model, leaving
the input model unchanged.

Inputs and Outputs
------------------
The primary driver of the master background step is a "spec3" type Association (ASN) file
or a ``ModelContainer`` populated from an ASN file. This is the ASN file used as input to
the :ref:`calwebb_spec3 <calwebb_spec3>` pipeline, which defines a stage 3 combined product
and its input members. The list of input members includes both "science" and "background"
exposure types. The master background subtraction step uses the input members designated
with "exptype: background" as the data from which to create the master background spectrum.
These should be :ref:`x1d <x1d>` products created for individual exposures at the end of
the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.
The master background will be subtracted from all input members designated as
"exptype: science" in the ASN, resulting in a new version of each such input. These inputs
should be :ref:`cal <cal>` products created for individual exposures by the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

For observing scenarios that use exposures of dedicated background targets, the "background"
members in the ASN will be the exposures from those background targets. For scenarios that
use source nodding within the field to sample the background, the lists of "science" and
"background" members in the ASN will use the same exposures. The only difference is that
the "science" members will refer to :ref:`cal <cal>` products, while the "background"
members will refer to :ref:`x1d <x1d>` products.

For example, the ASN file for a simple 2-point nodded observation consisting of two
exposures looks like the following::

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

As you can see from the above ASN list, the same two exposures are defined as
being both "science" and "background" members, because they both contain the target
of interest and a region of background. The "science" members, which are the
:ref:`cal <cal>` products created by the :ref:`calwebb_spec2 <calwebb_spec2>`
pipeline, are the data files that will have the master background subtraction
applied, while the "background" members are :ref:`x1d <x1d>` 1-D spectral
products from which the master background spectrum will be created.

Creating the 1-D master background spectrum
-------------------------------------------
The 1-D master background spectrum is created by combining data contained in the
:ref:`x1d <x1d>` products listed in the input ASN as being "exptype: background" members.
As noted above, the background members can be exposures obtained of dedicated background targets
or can be a collection of exposures of a point-like source observed in a nod pattern
(e.g. MIRI LRS fixed-slit "ALONG-SLIT-NOD" or NIRSpec IFU "2-POINT-NOD" dither patterns).

For the case of dedicated background target exposures, the 1-D spectrum contained in the
"FLUX" column of the background :ref:`x1d <x1d>` products will be used for creating the
master background spectrum. For the case of nodded exposures, the 1-D spectrum contained
in the "BACKGROUND" column of the :ref:`x1d <x1d>` products will be used.

When all the input background spectra have been collected, they are combined using the
:ref:`combine_1d <combine_1d_step>` step to produce the 1-D master background spectrum.
Because each input spectrum was originally created as the sum over a number of pixels
at a given wavelength, each spectrum is properly rescaled to yield background per
pixel before being combined.

Subtracting the master background
---------------------------------
The 1-D master background spectrum is interpolated by wavelength at each pixel of a 2-D source
spectrum and subtracted from it. The source data instances can be, for example, a set
of NIRSpec or MIRI IFU exposures, a set of NIRSpec MOS or fixed-slit 2-D extractions, or a set of
nodded MIRI LRS fixed-slit exposures. The subtraction process performs a loop over all input
source data instances and for each one it does the following:

 - Compute a 2-D wavelength grid corresponding to the 2-D source data. For some observing modes,
   such as NIRSpec MOS and fixed-slit, a 2-D wavelength array is computed and attached to the data
   in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline :ref:`extract_2d <extract_2d_step>` step.
   If such a wavelength array is present, it is used. For modes that don't have a 2-D
   wavelength array contained in the data product, it is computed on the fly using the WCS object
   for each source data instance.

 - Compute the background signal at each pixel in the 2-D wavelength grid by interpolating within
   the 1-D master background spectrum as a function of wavelength.

 - Subtract the resulting 2-D background image from the 2-D source data.

