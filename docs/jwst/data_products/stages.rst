Processing Levels and Product Stages
====================================
Here we describe the structure and content of the most frequently used forms of files for
JWST science data products, the vast majority of which are in FITS format. Each type of FITS
file is the result of serialization of a corresponding data model. All
JWST pipeline input and output products, with the exception of a few
reference files and catalogs, are serialized as FITS files.
The `ASDF <https://asdf-standard.readthedocs.io/en/stable/>`_ representation
of the data model is serialized as a FITS BINTABLE extension
within the FITS file, with EXTNAME="ASDF". The ASDF extension is essentially a
text character serialization in `YAML <https://yaml.org>`_ format of the
data model. The ASDF representation is read from the extension when a FITS file
is loaded into a data model to provide the initial instance of the data model.
Values in the other FITS extensions then either override this initial model or are added to it.

Within the various STScI internal data processing and archiving systems that are used for
routine processing of JWST data, there are some different uses of terminology to refer to
different levels or stages of processing and products. For those who are interested or
need to know, the table below gives high-level translations between those naming conventions.

+----------------------------------+-------------------------------------+------------------------------------+
| Data Processing Levels           | User Data Product Stages            | MAST/CAOM Data Levels              |
+==================================+=====================================+====================================+
| N/A                              | N/A                                 | -1 = Planned, but not yet executed |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 0 = Science telemetry      | Not available to users              | Not available to users             |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 0.5 = POD files            | Not available to users              | Not available to users             |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 1a = Original FITS file    | Stage 0 = Original FITS file        | 0 = raw                            |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 1b = Uncal FITS file       | Stage 0 = Fully-populated FITS file | 1 = uncalibrated                   |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 2a = Countrate exposure    | Stage 1 = Countrate FITS file       | 2 = calibrated                     |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 2b = Calibrated exposure   | Stage 2 = Calibrated exposure       | 2 = calibrated                     |
|                                  |                                     |                                    |
| \*Level 2c = CR-flagged exposure |                                     |                                    |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 3 = Combined data          | Stage 3 = Combined data             | 3 = Science product                |
+----------------------------------+-------------------------------------+------------------------------------+
| Level 4 = Analysis results       | Stage 4 = High-level product        | 4 = Contributed product            |
+----------------------------------+-------------------------------------+------------------------------------+

\*Note that Level 2c files are intermediate files produced during pipeline Stage 3 processing,
and are not final products (as opposed to all the other product types that are listed here).
Therefore, Level 2c files are not a final product of any pipeline stage, but are produced
within the pipeline Stage 3 processing. Level 2c files (identified by the 'crf' extension)
are in the same format as Level 2b products, with the difference being that their data quality
flags have been updated after running outlier detection in pipeline Stage 3 processing.

Throughout this document, we will use the "Stage" terminology to refer to data products.
Stage 0, 1, and 2 products are always files containing the data from a single exposure and a
single detector. A NIRCam exposure that uses all 10 detectors will therefore result in 10 separate
FITS files for the Stage 0, 1, and 2 products. Because these stages contain the data for a single
exposure, they are refered to as "exposure-based" products and use an "exposure-based" file naming
syntax. Stage 3 and 4 products, on the other hand, are constructed from the combined data of
multiple exposures for a given source or target. They are referred to as "source-based" products
and use a "source-based" file naming syntax. Observing modes that include multiple defined sources
within a single exposure or observation, such as NIRSpec MOS and NIRCam/NIRISS WFSS, will result in
multiple Stage 3 products, one for each defined or identifiable source.

