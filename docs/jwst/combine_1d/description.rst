Description
===========
The combine_1d step computes a weighted average of 1-D spectra and writes
the combined 1-D spectrum as output.

The weight will typically be the integration time or the exposure time,
but uniform (unit) weighting can be specified instead.

The combination of spectra proceeds as follows.  For each pixel of each
input spectrum, the corresponding pixel in the output is identified
(based on wavelength), and the input value multiplied by the weight is
added to the output buffer.  After all input spectra have been included,
the output is normalized by dividing by the sum of the weights.

The only part of this step that is not completely straightforward is the
determination of wavelengths for the output spectrum.  The output
wavelengths will be increasing, regardless of the order of the input
wavelengths.  In the ideal case, all input spectra would have wavelength
arrays that were very nearly the same.  In this case, each output
wavelength would be computed as the average of the wavelengths at the same
pixel in all the input files.  The combine_1d step is intended to handle a
more general case where the input wavelength arrays may be offset with
respect to each other, or they might not align well due to different
distortions.

All the input wavelength arrays will be concatenated and then sorted.
The code then looks for "clumps" in wavelength, based on the standard
deviation of a slice of the concatenated and sorted array of input
wavelengths; a small standard deviation implies a clump.  In regions of
the spectrum where the input wavelengths overlap with somewhat random
offsets and don't form any clumps, the output wavelengths are computed
as averages of the concatenated, sorted input wavelengths taken N at a
time, where N is the number of overlapping input spectra at that point.

Input
=====
An association file specifies which file or files to read for the input
data.  Each input data file contains one or more 1-D spectra in table
format, e.g. as written by the extract_1d step.  An input data file can
be either SpecModel (for one spectrum) or MultiSpecModel format (which
can contain more than one spectrum).

The association file should have an object called "products", which is
a one-element list containing a dictionary.  This dictionary contains two
entries (at least), one with key "name" and one with key "members".  The
value for key "name" is a string, the name that will be used for the output
file.  "members" is a list of dictionaries.  Each of these dictionaries
contains one input file name, identified by key "expname".

Output
======
The output will be in CombinedSpecModel format, with a table extension
having the name COMBINE1D.  This extension will have columns giving the
wavelength, countrate in counts/second, the sum of the weights that were
used when combining the input spectra, and the number of input spectra
that contributed to each output pixel.
