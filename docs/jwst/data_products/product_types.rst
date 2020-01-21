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

+--------------------+-----------------------+----------------------------+-------+--------+-------------------+---------------------------------------+
| Pipeline           | Input                 |  Output                    | Stage | base   | Units             | Description                           |
+====================+=======================+============================+=======+========+===================+=======================================+
| N/A                |                       | :ref:`uncal <uncal>`       |   0   | Exp    | DN                | Uncalibrated 4-D exposure             |
+--------------------+-----------------------+----------------------------+-------+--------+-------------------+---------------------------------------+
| calwebb_dark       | :ref:`uncal <uncal>`  | :ref:`dark <dark>`         |   1   | Exp    | DN                | 4-D corrected dark exposure           |
+--------------------+-----------------------+----------------------------+-------+--------+-------------------+---------------------------------------+
| calwebb_detector1  | :ref:`uncal <uncal>`  | :ref:`trapsfilled <trfld>` |   1   | Exp    | N/A               | Charge trap state data                |
|                    |                       +----------------------------+       |        +-------------------+---------------------------------------+
|                    |                       | :ref:`rateints <rate>`     |       |        | DN/s              | 3-D countrate data, per integration   |
|                    |                       +----------------------------+       |        |                   +---------------------------------------+
|                    |                       | :ref:`rate <rate>`         |       |        |                   | 2-D countrate data, per exposure      |
|                    |                       +----------------------------+       |        +-------------------+---------------------------------------+
|                    |                       | :ref:`fitopt <fitopt>`     |       |        | various           | *Optional* products from ramp_fit step|
|                    |                       +----------------------------+       |        +-------------------+---------------------------------------+
|                    |                       | :ref:`dark <dark>`         |       |        | DN                | *Optional* 4-D corrected dark exposure|
|                    |                       +----------------------------+       |        |                   +---------------------------------------+
|                    |                       | :ref:`ramp <ramp>`         |       |        |                   | *Optional* 4-D corrected ramp data    |
+--------------------+-----------------------+----------------------------+-------+--------+-------------------+---------------------------------------+

Stage 2 Data Products
+++++++++++++++++++++

+--------------------+-----------------------+----------------------------+-------+-------------------+-----------------------------------------+
| Pipeline           | Input                 |  Output                    | Base  | Units             | Description                             |
+====================+=======================+============================+=======+===================+=========================================+
| calwebb_image2     | :ref:`rateints <rate>`| :ref:`calints <cal>`       | Exp   | MJy/sr, MJy [#1]_ | | 3-D calibrated data,                  |
|                    |                       |                            |       |                   | | coronagraphys or TSO                  |
|                    +-----------------------+----------------------------|       |                   +-----------------------------------------+
| calwebb_tso-image2 | :ref:`rate <rate>`    | :ref:`cal <cal>`           |       |                   | | 2-D calibrated data                   |
|                    |                       +----------------------------|       |                   +-----------------------------------------+
|                    |                       | :ref:`bsub <bsub>`         |       |                   | | 2-D background-subtracted data,       |
|                    |                       |                            |       |                   | | when step applied                     |
|                    |                       +----------------------------|       |                   +-----------------------------------------+
|                    |                       | :ref:`i2d <i2d>`           |       |                   | | 2-D resampled imaging data            |
+--------------------+-----------------------+----------------------------+-------+-------------------+-----------------------------------------+
| calwebb_spec2      | :ref:`rateints <rate>`| :ref:`calints <cal>`       | Exp   | MJy/sr, MJy [#1]_ | | 3-D calibrated data,                  |
|                    |                       |                            |       |                   | | coronagraphys or TSO                  |
|                    |                       +----------------------------|       +-------------------+-----------------------------------------+
| calwebb_tso-spec2  |                       | :ref:`x1dints <x1d>`       |       | various           | | 1-D spectral data, per integration    |
|                    |                       +----------------------------|       +-------------------+-----------------------------------------+
|                    |                       | :ref:`flat <flat>`         |       | N/A               | | *Optional* for NIRSpec data.          |
|                    |                       |                            |       |                   | | Extracted combined flat.              |
|                    +-----------------------+----------------------------+-------+-------------------+-----------------------------------------+
|                    | :ref:`rate <rate>`    | :ref:`cal <cal>`           | Exp   | MJy/sr, MJy [#1]_ | | 2-D calibrated data                   |
|                    |                       +----------------------------|       |                   +-----------------------------------------+
|                    |                       | :ref:`bsub <bsub>`         |       |                   | | 2-D background-subtracted data        |
|                    |                       |                            |       |                   | | when step applied [#2]_               |
|                    |                       +----------------------------|       |                   +-----------------------------------------+
|                    |                       | :ref:`s3d <s3d>`           |       |                   | | 3-D resampled spectroscopic data.     |
|                    |                       |                            |       |                   | | For NRS_IFU or MIR_MRS                |
|                    |                       +----------------------------|       |                   +-----------------------------------------+
|                    |                       | :ref:`s2d <s2d>`           |       |                   | | 2-D resampled spectroscopic data      |
|                    |                       +----------------------------|       +-------------------+-----------------------------------------+
|                    |                       | :ref:`x1d <x1d>`           |       | various           | | 1-D extracted spectral data           |
|                    |                       +----------------------------|       +-------------------+-----------------------------------------+
|                    |                       | :ref:`flat <flat>`         |       | N/A               | | *Optional* for NIRSpec data.          |
|                    |                       |                            |       |                   | | Extracted combined flat.              |
+--------------------+-----------------------+----------------------------+-------+-------------------+-----------------------------------------+

Stage 3 Data Products
+++++++++++++++++++++

+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| Pipeline           | Input                 |  Output                    | Base | Units             | | Description                           |
+====================+=======================+============================+======+===================+=========================================+
| calwebb_image3     | :ref:`rateints <rate>`| :ref:`crfints <crf>`       | Exp  | MJy/sr, MJy [#1]_ | | 3-D CR-flagged calibrated data,       |
|                    |                       |                            |      |                   | | cornagraphy, TSO                      |
|                    +-----------------------+----------------------------+      |                   +-----------------------------------------+
|                    | :ref:`rate <rate>`    | :ref:`crf <crf>`           |      |                   | | 2-D CR-flagged calibrated data        |
|                    +-----------------------+----------------------------+------+                   +-----------------------------------------+
|                    | :ref:`cal <cal>`      | :ref:`i2d <i2d>`           | Src  |                   | | 2-D resampled imaging data            |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`cat <cat>`           |      |                   | | Source catalog                        |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| calwebb_spec3      | :ref:`rate <rate>`    | :ref:`crf <crf>`           | Exp  | MJy/sr, MJy [#1]_ | | 2-D CR-flagged calibrated data        |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    | :ref:`rateints <rate>`| :ref:`crfints <crf>`       |      |                   | | 3-D CR-flagged calibrated data,       |
|                    |                       |                            |      |                   | | cornagraphy, TSO                      |
|                    +-----------------------+----------------------------+------+                   +-----------------------------------------+
|                    | :ref:`cal <cal>`      | :ref:`s2d <s2d>`           | Src  |                   | | 2-D resampled spectroscopic data.     |
|                    |                       |                            |      |                   | | No IFU                                |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`s3d <s3d>`           |      |                   | | 3-D resampled spectroscopic data.     |
|                    |                       |                            |      |                   | | For NRS_IFU or MIR_MRS                |
|                    |                       +----------------------------+      +-------------------+-----------------------------------------+
|                    |                       | :ref:`x1d <x1d>`           |      | various           | | 1-D extracted spectroscopic data      |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| calwebb_ami3       | :ref:`cal <cal>`      | :ref:`ami <ami>`           | Src  | MJy/sr, MJy [#1]_ | | Fringe parameters                     |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`amiavg <ami>`        |      |                   | | Averaged fringe parameters            |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`aminorm <ami>`       |      |                   | | Normalized fringe parameters          |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| calwebb_coron3     | :ref:`rateints <rate>`| :ref:`crfints <crf>`       | Exp  | MJy/sr, MJy [#1]_ | | 3-D CR-flagged calibrated data        |
|                    +-----------------------+----------------------------+------+                   +-----------------------------------------+
|                    | :ref:`calints <cal>`  | :ref:`psfstack <psfstack>` | Src  |                   | | PSF library images                    |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`psfalign <psfalign>` | Exp  |                   | | Aligned PSF images                    |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`psfsub <psfsub>`     | Exp  |                   | | PSF-subtracted images                 |
|                    |                       +----------------------------+      |                   +-----------------------------------------+
|                    |                       | :ref:`i2d <i2d>`           | Src  |                   | | 2-D resampled PSF-subtracted image    |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| calwebb_tso3       | :ref:`rateints <rate>`| :ref:`crfints <crfints>`   | Exp  | MJy/sr, MJy [#1]_ | | 3-D CR-flagged calibrated data        |
|                    +-----------------------+----------------------------+------+-------------------+-----------------------------------------+
|                    | :ref:`calints <cal>`  | :ref:`phot <phot>`         | Src  | mag               | | TSO imaging photometry catalog        |
|                    |                       +----------------------------+      +-------------------+-----------------------------------------+
|                    |                       | :ref:`x1dints <x1dints>`   |      | various           | | TSO 1-D extracted spectra             |
|                    |                       +----------------------------+      +-------------------+-----------------------------------------+
|                    |                       | :ref:`whtlt <whtlt>`       |      | N/A               | | TSO spectral white-light catalog      |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+
| calwebb_wfs-image3 |  :ref:`cal <cal>`     | :ref:`wfscmb <wfscmb>`     | Src  | MJy/sr, MJy [#1]_ | | 2-D combined WFS&C image              |
+--------------------+-----------------------+----------------------------+------+-------------------+-----------------------------------------+

.. rubric :: Footnotes

.. [#1] NIRSpec and NIRISS SOSS point sources have MJy units
.. [#2] Non-TSO. To be implemented 
