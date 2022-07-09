===================
Available Pipelines
===================
There are many pre-defined pipeline modules for processing
data from different instrument observing modes through each of the 3 stages
of calibration. For all of the details see :ref:`pipelines`.

.. _pipeline_step_suffix_definitions:

Pipeline/Step Suffix Definitions
================================

However the output file name is determined (:ref:`see above
<intro_output_file_discussion>`), the various stage 1, 2, and 3 pipeline modules
will use that file name, along with a set of predetermined suffixes, to compose
output file names. The output file name suffix will always replace any known
suffix of the input file name. Each pipeline module uses the appropriate suffix
for the product(s) it is creating. The list of suffixes is shown in the
following table. Replacement occurs only if the suffix is one known to the
calibration code. Otherwise, the new suffix will simply be appended to the
basename of the file.

=============================================  ========
Product                                        Suffix
=============================================  ========
Uncalibrated raw input                         uncal
Corrected ramp data                            ramp
Corrected countrate image                      rate
Corrected countrate per integration            rateints
Optional fitting results from ramp_fit step    fitopt
Background-subtracted image                    bsub
Per integration background-subtracted image    bsubints
Calibrated image                               cal
Calibrated per integration images              calints
CR-flagged image                               crf
CR-flagged per integration images              crfints
Resampled 2D image                             i2d
Resampled 2D spectrum                          s2d
Resampled 3D IFU cube                          s3d
1D extracted spectrum                          x1d
1D extracted spectra per integration           x1dints
1D combined spectrum                           c1d
Source catalog                                 cat
Segmentation map                               segm
Time Series photometric catalog                phot
Time Series white-light catalog                whtlt
Coronagraphic PSF image stack                  psfstack
Coronagraphic PSF-aligned images               psfalign
Coronagraphic PSF-subtracted images            psfsub
AMI fringe and closure phases                  ami
AMI averaged fringe and closure phases         amiavg
AMI normalized fringe and closure phases       aminorm
=============================================  ========


For More Information
====================
More information on logging and running pipelines can be found in the ``stpipe``
User's Guide at :ref:`stpipe-user-steps`.

More detailed information on writing pipelines can be found
in the ``stpipe`` Developer's Guide at :ref:`stpipe-devel-steps`.

If you have questions or concerns regarding the software, please open an issue
at https://github.com/spacetelescope/jwst/issues or contact
the `JWST Help Desk <https://jwsthelp.stsci.edu>`_.