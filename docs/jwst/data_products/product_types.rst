Data Product Types
------------------
The following table contains a list of all data product types, as given by their file name suffix,
the pipeline and stage of processing that creates each type, and an indication of whether the file naming
is exposure-based or source-based. The product name suffixes are active links to detailed descriptions
in the following sections.

+--------------------+----------------------------+-------+------------+-------------------------------------+
| Pipeline           | Type Suffix                | Stage | Exp/Source | Description                         |
+====================+============================+=======+============+=====================================+
| N/A                | :ref:`uncal <uncal>`       |   0   | Exposure   | Uncalibrated 4-D exposure           |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_dark       | :ref:`dark <dark>`         |   1   | Exposure   | 4-D corrected dark exposure         |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_detector1  | :ref:`trapsfilled <trfld>` |   1   | Exposure   | Charge trap state data              |
|                    |                            |       |            |                                     |
|                    | :ref:`ramp <ramp>`         |   1   | Exposure   | 4-D corrected ramp data             |
|                    |                            |       |            |                                     |
|                    | :ref:`rateints <rate>`     |   1   | Exposure   | 3-D countrate data, per integration |
|                    |                            |       |            |                                     |
|                    | :ref:`rate <rate>`         |   1   | Exposure   | 2-D countrate data, per exposure    |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_spec2      | :ref:`bsubints <bsub>`     |   2   | Exposure   | 3-D background-subtracted data      |
|                    |                            |       |            |                                     |
|                    | :ref:`bsub <bsub>`         |   2   | Exposure   | 2-D background-subtracted data      |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_image2     | :ref:`calints <cal>`       |   2   | Exposure   | 3-D calibrated data                 |
|                    |                            |       |            |                                     |
| calwebb_tso-image2 | :ref:`cal <cal>`           |   2   | Exposure   | 2-D calibrated data                 |
|                    |                            |       |            |                                     |
| calwebb_wfs-image2 | :ref:`i2d <i2d>`           |   2   | Exposure   | 2-D resampled imaging data          |
|                    |                            |       |            |                                     |
| calwebb_spec2      | :ref:`s2d <s2d>`           |   2   | Exposure   | 2-D resampled spectroscopic data    |
|                    |                            |       |            |                                     |
| calwebb_tso-spec2  | :ref:`s3d <s3d>`           |   2   | Exposure   | 3-D resampled spectroscopic data    |
|                    |                            |       |            |                                     |
|                    | :ref:`x1dints <x1d>`       |   2   | Exposure   | 1-D spectral data, per integration  |
|                    |                            |       |            |                                     |
|                    | :ref:`x1d <x1d>`           |   2   | Exposure   | 1-D extracted spectral data         |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_image3     | :ref:`crf <crf>`           |   2   | Exposure   | 2-D CR-flagged calibrated data      |
|                    |                            |       |            |                                     |
|                    | :ref:`i2d <i2d>`           |   3   | Source     | 2-D resampled imaging data          |
|                    |                            |       |            |                                     |
|                    | :ref:`cat <cat>`           |   3   | Source     | Source catalog                      |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_spec3      | :ref:`crf <crf>`           |   2   | Exposure   | 2-D CR-flagged calibrated data      |
|                    |                            |       |            |                                     |
|                    | :ref:`s2d <s2d>`           |   3   | Source     | 2-D resampled spectroscopic data    |
|                    |                            |       |            |                                     |
|                    | :ref:`s3d <s3d>`           |   3   | Source     | 3-D resampled spectroscopic data    |
|                    |                            |       |            |                                     |
|                    | :ref:`x1d <x1d>`           |   3   | Source     | 1-D extracted spectroscopic data    |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_ami3       | :ref:`ami <ami>`           |   3   | Source     | Fringe parameters                   |
|                    |                            |       |            |                                     |
|                    | :ref:`amiavg <ami>`        |   3   | Source     | Averaged fringe parameters          |
|                    |                            |       |            |                                     |
|                    | :ref:`aminorm <ami>`       |   3   | Source     | Normalized fringe parameters        |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_coron3     | :ref:`crfints <crf>`       |   2   | Exposure   | 3-D CR-flagged calibrated data      |
|                    |                            |       |            |                                     |
|                    | :ref:`psfstack <psfstack>` |   3   | Source     | PSF library images                  |
|                    |                            |       |            |                                     |
|                    | :ref:`psfalign <psfalign>` |   3   | Exposure   | Aligned PSF images                  |
|                    |                            |       |            |                                     |
|                    | :ref:`psfsub <psfsub>`     |   3   | Exposure   | PSF-subtracted images               |
|                    |                            |       |            |                                     |
|                    | :ref:`i2d <i2d>`           |   3   | Source     | 2-D resampled PSF-subtracted image  |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_tso3       | :ref:`crfints <crfints>`   |   2   | Exposure   | 3-D CR-flagged calibrated data      |
|                    |                            |       |            |                                     |
|                    | :ref:`phot <phot>`         |   3   | Source     | TSO imaging photometry catalog      |
|                    |                            |       |            |                                     |
|                    | :ref:`x1dints <x1dints>`   |   3   | Source     | TSO 1-D extracted spectra           |
|                    |                            |       |            |                                     |
|                    | :ref:`whtlt <whtlt>`       |   3   | Source     | TSO spectral white-light catalog    |
+--------------------+----------------------------+-------+------------+-------------------------------------+
| calwebb_wfs-image3 | :ref:`wfscmb <wfscmb>`     |   3   | Source     | 2-D combined WFS&C image            |
+--------------------+----------------------------+-------+------------+-------------------------------------+

