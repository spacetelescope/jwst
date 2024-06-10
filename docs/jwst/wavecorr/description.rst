Description
============

:Class: `jwst.wavecorr.WavecorrStep`
:Alias: wavecorr

The wavelength correction (``wavecorr``) step in the
:ref:`calwebb_spec2 <calwebb_spec2>` pipeline updates the wavelength
assignments for NIRSpec fixed-slit (FS) and MOS point sources that are
known to be off center (in the dispersion direction) in their slit.

NIRSpec MOS
-----------
For NIRSpec MOS exposures (EXP_TYPE="NRS_MSASPEC"), wavelength
assignments created during :ref:`extract_2d <extract_2d_step>` are based on
a source that's perfectly centered in a slitlet. Most sources, however,
are not centered in every slitlet in a real observation.
The MSA meta data assigned to each slitlet in the
:ref:`extract_2d <extract_2d_step>` step includes estimates of the source
x (dispersion) and y (cross-dispersion) location within the slitlet.
These are recorded in the "SRCXPOS" and "SRCYPOS" keywords in the SCI
extension header of each slitlet in a FITS product.

The ``wavecorr`` step loops over all slit instances in the input
science product and updates the WCS models of slits that contain a point 
source to include a wavelength correction. The point source determination is 
based on the value of the "SRCTYPE" keyword populated for each slit by the
:ref:`srctype <srctype_step>` step. The computation of the correction is
based on the "SRCXPOS" value. A value of 0.0 indicates a perfectly centered
source, and ranges from -0.5 to +0.5 for sources at the extreme edges
of a slit. The computation uses calibration data from the ``WAVECORR``
reference file, which contains pixel shifts as a function of source position 
and wavelength, and can be converted to wavelength shifts with the dispersion. 
For each slit, the ``wavecorr`` step uses the average wavelengths and 
dispersions in a slit (averaged across the cross-dispersion direction) to 
calculate corresponding corrected wavelengths.  It then uses the average 
wavelengths and their corrections to generate a transform that interpolates 
between "center of slit" wavelengths and corrected wavelengths.  This 
transformation is added to the slit WCS after the ``slit_frame`` and
produces a new wavelength corrected slit frame, ``wavecorr_frame``.

NIRSpec Fixed Slit (FS)
-----------------------
Fixed slit data do not have an *a priori* estimate of the source
location within a given slit, so the estimated source location is
computed by the ``wavecorr`` step. It uses the target coordinates in
conjunction with the aperture reference point in V2/V3 space to
estimate the fractional location of the source within the given slit.
Note that this computation can only be performed for the primary slit
in the exposure, which is given in the "FXD_SLIT" keyword. The positions
of sources in any additional slits cannot be estimated and therefore
the wavelength correction is only applied to the primary slit.\ :sup:`1`

The estimated position of the source within the primary slit (in the
dispersion direction) is then used in the same manner as described above
for MOS slitlets to update the slit WCS and compute corrected wavelengths
for the primary slit.

Upon successful completion of the step, the status keyword "S_WAVCOR"
is set to "COMPLETE".

:sup:`1`\ Note that fixed slits that are planned as part of a combined
MOS and FS observation do have *a priori* estimates of their source
locations, via the :ref:`MSA metadata file<msa_metadata>`. When available,
these source locations are directly used, instead of recomputing the source
position in the ``wavecorr`` step.
