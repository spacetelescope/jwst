Data Product Types
------------------
The following tables contain lists of all data product types, as given by their file name suffix. There is one table per stage of processing.
All tables indicate whether the file naming is exposure-based (Exp) or source-based (Src).
When the product is not created by default, the flag *Optional* is indicated in the
description. The different stages of the calibration pipeline are as defined in
the `Algorithms Documentation <https://jwst-docs.stsci.edu/jwst-data-reduction-pipeline/algorithm-documentation>`_.
The product name suffixes are active links to detailed descriptions in the following sections.

Stage 0 and Stage 1 Data Products
+++++++++++++++++++++++++++++++++

+----------------------------------------------+-----------------------+----------------------------+-------+------+---------+----------------------------------------+
| Pipeline                                     | Input                 |  Output(s)                 | Stage | Base | Units   | Description                            |
+==============================================+=======================+============================+=======+======+=========+========================================+
| N/A                                          |                       | :ref:`uncal <uncal>`       |   0   | Exp  | DN      | Uncalibrated 4-D exposure data         |
+----------------------------------------------+-----------------------+----------------------------+-------+------+---------+----------------------------------------+
| :ref:`calwebb_dark <calwebb_dark>`           | :ref:`uncal <uncal>`  | :ref:`dark <dark>`         |   1   | Exp  | DN      | 4-D corrected dark exposure data       |
+----------------------------------------------+-----------------------+----------------------------+-------+------+---------+----------------------------------------+
| :ref:`calwebb_detector1 <calwebb_detector1>` | :ref:`uncal <uncal>`  | :ref:`trapsfilled <trfld>` |   1   | Exp  | N/A     | Charge trap state data                 |
|                                              |                       +----------------------------+       |      +---------+----------------------------------------+
|                                              |                       | :ref:`rateints <rateints>` |       |      | DN/s    | 3-D countrate data (per integration)   |
|                                              |                       +----------------------------+       |      |         +----------------------------------------+
|                                              |                       | :ref:`rate <rate>`         |       |      |         | 2-D countrate data (per exposure)      |
|                                              |                       +----------------------------+       |      +---------+----------------------------------------+
|                                              |                       | fitopt                     |       |      | various | *Optional* fit info from ramp_fit step |
|                                              |                       +----------------------------+       |      +---------+----------------------------------------+
|                                              |                       | dark                       |       |      | DN      | *Optional* 3-D on-the-fly dark data    |
|                                              |                       +----------------------------+       |      |         +----------------------------------------+
|                                              |                       | :ref:`ramp <ramp>`         |       |      |         | *Optional* 4-D corrected ramp data     |
+----------------------------------------------+-----------------------+----------------------------+-------+------+---------+----------------------------------------+

Stage 2 Data Products
+++++++++++++++++++++

+----------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+
| Pipeline                               | Input                  |  Output(s)               | Base | Units                 | Description                           |
+========================================+========================+==========================+======+=======================+=======================================+
| :ref:`calwebb_image2 <calwebb_image2>` | :ref:`rate <rate>`     | :ref:`bsub <bsub>`       | Exp  | DN/s                  | | 2-D background-subtracted data,     |
|                                        |                        |                          |      |                       | | when background step applied        |
|                                        |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | :ref:`cal <cal>`         |      | MJy/sr, MJy [#1]_     | | 2-D calibrated data                 |
|                                        |                        +--------------------------+      |                       +---------------------------------------+
|                                        |                        | :ref:`i2d <i2d>`         |      |                       | | 2-D resampled imaging data          |
+----------------------------------------+------------------------+--------------------------+      +-----------------------+---------------------------------------+
| :ref:`calwebb_image2 <calwebb_image2>` | :ref:`rateints <rate>` | :ref:`calints <calints>` |      | MJy/sr, MJy [#1]_     | | 3-D calibrated data;                |
| with TSO data                          |                        |                          |      |                       | | coronagraphy and TSO                |
+----------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+
| :ref:`calwebb_spec2 <calwebb_spec2>`   | :ref:`rate <rate>`     | :ref:`bsub <bsub>`       | Exp  | DN/s                  | | 2-D background-subtracted data,     |
|                                        |                        |                          |      |                       | | when background step applied        |
|                                        |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | :ref:`cal <cal>`         |      | MJy/sr, MJy [#1]_     | | 2-D calibrated data                 |
|                                        |                        +--------------------------+      |                       +---------------------------------------+
|                                        |                        | :ref:`s3d <s3d>`         |      |                       | | 3-D resampled spectroscopic data;   |
|                                        |                        |                          |      |                       | | NIRSpec IFU and MIRI MRS            |
|                                        |                        +--------------------------+      |                       +---------------------------------------+
|                                        |                        | :ref:`s2d <s2d>`         |      |                       | | 2-D resampled spectroscopic data    |
|                                        |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | :ref:`x1d <x1d>`         |      | various               | | 1-D extracted spectral data         |
|                                        |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | flat                     |      | N/A                   | | *Optional* for NIRSpec data;        |
|                                        |                        |                          |      |                       | | on-the-fly constructed flat.        |
+----------------------------------------+------------------------+--------------------------+      +-----------------------+---------------------------------------+
| :ref:`calwebb_spec2 <calwebb_spec2>`   | :ref:`rateints <rate>` | :ref:`calints <calints>` |      | MJy/sr, MJy [#1]_     | | 3-D calibrated data; TSO            |
| with TSO data                          |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | :ref:`x1dints <x1dints>` |      | various               | | 1-D spectral data (per integration) |
|                                        |                        +--------------------------+      +-----------------------+---------------------------------------+
|                                        |                        | flat                     |      | N/A                   | | *Optional* for NIRSpec data;        |
|                                        |                        |                          |      |                       | | on-the-fly constructed flat.        |
+----------------------------------------+------------------------+--------------------------+------+-----------------------+---------------------------------------+

Stage 3 Data Products
+++++++++++++++++++++

+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| Pipeline                                       | Input                |  Outputs                   | Base | Units                 | | Description                        |
+================================================+======================+============================+======+=======================+======================================+
| :ref:`calwebb_image3 <calwebb_image3>`         | :ref:`cal <cal>`     | :ref:`crf <crf>`           | Exp  | MJy/sr, MJy [#1]_     | | 2-D CR-flagged calibrated data     |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`i2d <i2d>`           | Src  |                       | | 2-D resampled imaging data         |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`cat <cat>`           |      | N/A                   | | Source catalog                     |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`segm <segm>`         |      | N/A                   | | Segmentation map                   |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| :ref:`calwebb_spec3 <calwebb_spec3>`           | :ref:`cal <cal>`     | :ref:`crf <crf>`           | Exp  | MJy/sr, MJy [#1]_     | | 2-D CR-flagged calibrated data     |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`s2d <s2d>`           | Src  |                       | | 2-D resampled spectroscopic data;  |
|                                                |                      |                            |      |                       | | Non-IFU                            |
|                                                |                      +----------------------------+      |                       +--------------------------------------+
|                                                |                      | :ref:`s3d <s3d>`           |      |                       | | 3-D resampled spectroscopic data;  |
|                                                |                      |                            |      |                       | | NIRSpec IFU and MIRI MRS           |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`x1d <x1d>`           |      | various               | | 1-D extracted spectroscopic data   |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`c1d <c1d>`           |      | various               | | 1-D combined spectroscopic data    |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| :ref:`calwebb_ami3 <calwebb_ami3>`             | :ref:`cal <cal>`     | :ref:`ami <ami>`           | Exp  | N/A                   | | Fringe parameters (per exposure)   |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`amiavg <ami>`        | Src  |                       | | Averaged fringe parameters         |
|                                                |                      +----------------------------+      |                       +--------------------------------------+
|                                                |                      | :ref:`aminorm <ami>`       |      |                       | | Normalized fringe parameters       |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| :ref:`calwebb_coron3 <calwebb_coron3>`         | :ref:`calints <cal>` | :ref:`crfints <crf>`       | Exp  | MJy/sr, MJy [#1]_     | | 3-D CR-flagged calibrated data     |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`psfstack <psfstack>` | Src  |                       | | PSF library images                 |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`psfalign <psfalign>` | Exp  |                       | | Aligned PSF images                 |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`psfsub <psfsub>`     | Exp  |                       | | PSF-subtracted images              |
|                                                |                      +----------------------------+------+                       +--------------------------------------+
|                                                |                      | :ref:`i2d <i2d>`           | Src  |                       | | 2-D resampled PSF-subtracted image |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| :ref:`calwebb_tso3 <calwebb_tso3>`             | :ref:`calints <cal>` | :ref:`crfints <crfints>`   | Exp  | MJy/sr, MJy [#1]_     | | 3-D CR-flagged calibrated data     |
|                                                |                      +----------------------------+------+-----------------------+--------------------------------------+
|                                                |                      | :ref:`phot <phot>`         | Src  | N/A                   | | TSO imaging photometry catalog     |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`x1dints <x1dints>`   |      | various               | | TSO 1-D extracted spectra          |
|                                                |                      +----------------------------+      +-----------------------+--------------------------------------+
|                                                |                      | :ref:`whtlt <whtlt>`       |      | N/A                   | | TSO spectral white-light catalog   |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+
| :ref:`calwebb_wfs-image3 <calwebb_wfs-image3>` |  :ref:`cal <cal>`    | :ref:`wfscmb <wfscmb>`     | Src  | MJy/sr, MJy [#1]_     | | 2-D combined WFS&C image           |
+------------------------------------------------+----------------------+----------------------------+------+-----------------------+--------------------------------------+

.. [#1] NIRSpec and NIRISS SOSS point sources have MJy units; all others are MJy/sr
