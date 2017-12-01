
Description
===========
The background subtraction step in the calibration pipeline performs
image-from-image subtraction in order to accomplish subtraction of background
signal. The step takes as input one target exposure, to which the
subtraction will be applied, and a list of one or more background exposures.
Wide-field slitless spectroscopy data are processed differently from
all other types of data.  This is described below, toward the end of this
section.  Note, however, that background subtraction for WFSS data is
currently disabled, because there are no background reference files in CRDS.

If more than one background exposure is provided, they will be averaged
together before being applied to the target exposure.
The average background image is produced as follows:

 - The SCI arrays of all background exposures are averaged
 - The ERR arrays of all background exposures are summed in quadrature and
   then converted to an uncertainty in the mean
 - The DQ arrays of all background exposures are combined using a bitwise-OR
   operation

The average background exposure is then subtracted from the target exposure.
The subtraction consists of the following operations:

 - The SCI array of the average background is subtracted from the SCI
   array of the target

 - The ERR array of the target is not operated upon until full error
   propagation is implemented in the entire pipeline

 - The DQ arrays of the average background and the target are combined
   using a bitwise-OR operation


If the target exposure is a simple ImageModel, the background image is
subtracted from it. If the target exposure is in the form of a CubeModel
(e.g. the result of a time-series exposure), the background image
is subtracted from all planes of the CubeModel.

For wide-field slitless spectroscopy data (NIS_WFSS or NRC_GRISM), the
background reference image will be scaled to match the science data within
background regions, and the scaled image will be subtracted from the
science data.  The background regions in the science data are places where
there are no source spectra; the locations of the spectra are determined
from a source catalog (specified by primary header keyword SCATFILE)
and a reference file that gives the wavelength range (based on filter and
grating) that is relevant to the science data.  Robust mean values are
obtained for the background region in the science data and for the same
region in the background reference image, and the ratio of those two mean
values is used to scale the background reference image to the science data,
and the scaled background is subtracted.  The robust mean is computed by
excluding the lowest 25% and highest 25% of the data (using the
numpy.percentile function), and taking a simple arithmetic mean of the
remaining values.  Note that NaN values (if any) in the background
reference image are currently set to zero.  If there are a lot of NaNs,
it may be that more than 25% of the lowest values will need to be excluded.

The output results are always returned in a new data model, leaving the original
input model unchanged.

Upon successful completion of the step, the S_BKDSUB keyword will be set to
'COMPLETE' in the output product header.

Step Arguments
==============
There are no step-specific arguments.

Reference File
==============
The only mode for which this step uses reference files is wide-field
slitless spectroscopy.  For WFSS, there are two reference files.
Reference type "wfssbkg" is the actual background image, and reference
type "wavelengthrange" contains information about the range of wavelengths
in the exposure.  The latter is used, together with a source catalog, to
create a mask showing the locations of source spectra in the image, and
hence, where the background regions are.
