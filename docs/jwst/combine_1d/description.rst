Description
===========

:Class: `jwst.combine_1d.Combine1dStep`
:Alias: combine_1d

The ``combine_1d`` step computes a weighted average of 1-D spectra and writes
the combined 1-D spectrum as output.

The combination of spectra proceeds as follows.  For each pixel of each
input spectrum, the corresponding pixel in the output is identified
(based on wavelength), and the input value multiplied by the weight is
added to the output buffer.  Pixels that are flagged (via the DQ column)
with "DO_NOT_USE" will not contribute to the output.  After all input
spectra have been included, the output is normalized by dividing by
the sum of the weights.

The weight will typically be the integration time or the exposure time,
but uniform (unit) weighting can be specified instead.

The only part of this step that is not completely straightforward is the
determination of wavelengths for the output spectrum.  The output
wavelengths will be increasing, regardless of the order of the input
wavelengths.  In the ideal case, all input spectra would have wavelength
arrays that were very nearly the same.  In this case, each output
wavelength would be computed as the average of the wavelengths at the same
pixel in all the input files.  The combine_1d step is intended to handle a
more general case where the input wavelength arrays may be offset with
respect to each other, or they might not align well due to different
distortions.  All the input wavelength arrays will be concatenated and then
sorted.  The code then looks for "clumps" in wavelength, based on the
standard deviation of a slice of the concatenated and sorted array of input
wavelengths; a small standard deviation implies a clump.  In regions of
the spectrum where the input wavelengths overlap with somewhat random
offsets and don't form any clumps, the output wavelengths are computed
as averages of the concatenated, sorted input wavelengths taken N at a
time, where N is the number of overlapping input spectra at that point.

Input
=====
An association file specifies which file or files to read for the input
data.  Each input data file contains one or more 1-D spectra in table
format, e.g. as written by the extract_1d step.  Each input data file will
ordinarily be in MultiSpecModel format (which can contain more than one
spectrum).

The association file should have an object called "products", which is
a one-element list containing a dictionary.  This dictionary contains two
entries (at least), one with key "name" and one with key "members".  The
value for key "name" is a string, the name that will be used as a basis for
creating the output file name.  "members" is a list of dictionaries, each
of which contains one input file name, identified by key "expname".

Output
======
For most modes, the output will be in CombinedSpecModel format, with a table extension
having the name COMBINE1D.  This extension will have eight columns, giving
the wavelength, flux, error estimate for the flux, surface brightness,
error estimate for the surface brightness, the combined data quality flags,
the sum of the weights that were used when combining the input spectra,
and the number of input spectra that contributed to each output pixel.

For WFSS modes, which may have hundreds or thousands of spectra from different sources,
the output will be in WFSSMultiCombinedSpecModel format.
This model differs from the other MultiCombinedSpecModel classes in that
it is designed to hold all the spectra in a WFSS observation in a single
"flat" table format. Therefore, there is only one item per spectral order
in the `spec` list, and each object in the `spec` list has
a `spec_table` attribute that contains the spectral data and metadata
for all sources in the observation.

The spectral table for this model contains the same columns as the ``CombinedSpecModel``, but
each row in the table contains the combined spectrum for a single source. The spectral columns
are 2D: each row is a 1D vector containing all data points for the spectrum. In addition, the
spectral tables for this model have extra 1D columns to contain the metadata for the spectrum in each row.
These metadata fields include:
SOURCE_ID, N_ALONGDISP, SOURCE_TYPE, SOURCE_RA, SOURCE_DEC.

Note that the vector columns have the same length for all the sources in the table, meaning that
the number of elements in the table rows is set by the spectrum with the most data points.
The other spectra are NaN-padded to match the longest spectrum,
and the number of valid data points for each spectrum is recorded in the N_ALONGDISP column.

For example, to access the wavelength and flux for a specific source ID (say, 1200)
in a WFSSMultiCombinedSpecModel:

.. doctest-skip::

  >>> from stdatamodels.jwst import datamodels
  >>> model = datamodels.open('multi_wfss_c1d.fits')
  >>> spec_this_order = model.spec[0]
  >>> print(spec.spectral_order) # returns e.g. '1'
  >>> tab = spec.spec_table
  >>> row_want = tab[tab["SOURCE_ID"] == 1200][0]
  >>> nelem = row_want["N_ALONGDISP"]
  >>> wave, flux = row_want["WAVELENGTH"][:nelem], row_want["FLUX"][:nelem]
