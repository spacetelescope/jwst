Examples of How to Run Cube_Build ================================= It is
assumed that the input data have been processed through the
:ref:`calwebb_detector1 <calwebb_detector1>` pipeline and up through the
``photom`` step of the :ref:`calwebb_spec2 <calwebb_spec2>` pipeline.

Cube Building for MIRI data
-------------------------------
To run cube_build on a single MIRI exposure (containing channel 1 and 2), but only creating an IFU cube for channel 1::

  strun jwst.cube_build.CubeBuildStep MIRM103-Q0-SHORT_495_cal.fits --ch=1

The output 3D spectral cube will be saved in a file called MIRM103-Q0-SHORT_495_ch1-short_s3d.fits

To run cube_build using an association table containing 4 dithered images::

  strun jwst.cube_build.CubeBuildStep cube_build_4dither_asn.json

where the ASN file cube_build_4dither_asn.json contains::

	{"asn_rule": "Asn_MIRIFU_Dither",
         "target": "MYTarget",
         "asn_id": "c3001",
	 "asn_pool": "jw00024_001_01_pool",
         "program": "00024","asn_type":"dither",
	 "products": [
                     {"name": "MIRM103-Q0-Q3",
                     "members":
                      [{"exptype": "SCIENCE", "expname": "MIRM103-Q0-SHORT_495_cal.fits"},
                       {"exptype": "SCIENCE", "expname": "MIRM103-Q1-SHORT_495_cal.fits"},
                       {"exptype": "SCIENCE", "expname": "MIRM103-Q2-SHORT_495_cal.fits"},
                       {"exptype": "SCIENCE", "expname": "MIRM103-Q3-SHORT_495_cal.fits"}]}
	              ]
        }

The default output will be two IFU cubes. The first will contain the combined dithered images for
channel 1, sub-channel SHORT and the second will contain the channel 2, sub-channel SHORT data.
The output root file names are defined by the product "name" attribute in
the association table and results in files MIRM103-Q0-Q3_ch1-short_s3d.fits and MIRM103-Q0-Q3_ch2-short_s3d.fits.

To use the same association table, but combine all the data, use the output_type=multi option::

  strun jwst.cube_build.CubeBuildStep cube_build_4dither_asn.json --output_type=multi

The output IFU cube file will be MIRM103-Q0-Q3_ch1-2-short_s3d.fits


Cube building for NIRSpec data
----------------------------------

To run ``cube_build`` on a single NIRSpec exposure that uses grating G140H and filter F100LP::

  strun jwst.cube_build.CubeBuildStep jwtest1004001_01101_00001_nrs2_cal.fits

The output file will be jwtest1004001_01101_00001_nrs2_g140h-f100lp_s3d.fits

To run ``cube_build`` using an association table containing data from exposures using G140H+F100LP and G140H+F070LP::

  strun jwst.cube_build.CubeBuildStep nirspec_multi_asn.json

where the association file contains::

	{"asn_rule": "Asn_NIRSPECFU_Dither",
         "target": "MYTarget",
	 "asn_pool": "jw00024_001_01_pool",
	 "program": "00024","asn_type":"NRSIFU",
	 "asn_id":"a3001",
	 "products": [
         {"name": "JW3-6-NIRSPEC",
         "members":
         [{"exptype": "SCIENCE", "expname": "jwtest1003001_01101_00001_nrs1_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1004001_01101_00001_nrs2_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1005001_01101_00001_nrs1_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1006001_01101_00001_nrs2_cal.fits"}]}
         ]
	 }

The output will be two IFU cubes, one for each grating+filter combination: JW3-6-NIRSPEC_g140h-f070lp_s3d.fits and
JW3-6-NIRSPEC_g140h-f100lp_s3d.fits.

