Description
===========
The master spectroscopic background subtraction step subtracts background signal from
2-D spectroscopic data using a 1-D master background spectrum. The 1-D master background
spectrum is computed from one or more input exposures, or can alternatively be supplied
by the user. The 1-D background spectrum - flux versus wavelength - is projected into the
2-D space of source data based on the wavelength of each pixel in the 2-D data. The resulting
2-D background signal is then subtracted directly from the 2-D source data.

Upon successful completion of the step, the S_MSBSUB keyword is set to 'COMPLETE' in the
output product. The background-subtracted results are returned as a new data model, leaving
the input model unchanged.

Creating the 1-D master background spectrum
-------------------------------------------
The 1-D master background spectrum is created by combining the data contained in the
:ref:`x1d <x1d>` products created by the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline. The
collection of exposure products to be used as background is designated in the input ASN file,
by setting 'exptype=background' for the appropriate product members. The background members
could be exposures obtained of dedicated background targets or could be a collection of
on-target exposures for a point-source observed in a nod pattern (e.g. MIRI LRS fixed-slit
2-point nods).

The flux versus wavelength data in the :ref:`x1d <x1d>` products for the designated background members
is combined using the :ref:`combine_1d <combine_1d_step>` step to produce the 1-D master
background spectrum.

[Need a brief description here of how the multiple spectra are resampled/combined in
wavelength space to form the final MSB spectrum, once that's been decided.]

Subtracting the master background
---------------------------------
The 1-D master background spectrum is projected into the 2-D space of each source data instance
contained in the inputs and subtracted from it. The source data instances could be, for example, a set
of NIRSpec or MIRI IFU exposures, a set of NIRSpec MOS or fixed-slit 2-D extractions, or a set of
nodded MIRI LRS fixed-slit exposures. The subtraction process performs a loop over all input
source data instances and for each one it does the following:

 - Compute a 2-D wavelength grid corresponding to the 2-D source data. For some observing modes,
   such as NIRSpec MOS and fixed-slit, a 2-D wavelength is computed and attached to the data during
   the :ref:`extract_2d <extract_2d_step>` step in the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.
   If such a 2-D wavelength array is present, it is used. For other modes that don't have a 2-D
   wavelength array contained in the data product, it is computed on the fly using the WCS object
   for each source data instance.

 - Compute the background signal at each pixel in the 2-D wavelength grid by interpolating within
   the 1-D master background spectrum as a function of wavelength.

 - Subtract the resulting 2-D background image from the 2-D source data.

