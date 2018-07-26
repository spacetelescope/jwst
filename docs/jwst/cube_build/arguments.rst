.. _arguments:

Step Arguments
==============
The input data to the cube program can be a
single exposure, a data model passed from another pipeline step,  or a list of exposures contained in an association table.
The default spectral cubes contain data from a single band. For MIRI a single band covers a single channel and single sub-channel,
and for NIRSPEC a single band covers a single grating and single filter. There are options that allow the user to select the
type of data to be used in constructing the spectral IFU cube.
The arguments controlling the  types of data to use  are:

* ``--channel [integer]``

The only valid values for # are 1,2,3,4 or ALL .
This argument is only valid for MIRI data. If the ``channel`` argument is given, then only data corresponding to that channel
will be used in constructing the cube.  If the user wants to construct a cube from more than one channel,
then all the values are contained in the string with a comma between each channel number. For example,
to create a cube with channel 1 and 2 the argument list is ``--channel='1, 2'``. If this value is not set then the default is
to make single channel-single sub-channel IFU Cubes.

* ``--band [string]``

This is a MIRI option and the  only valid values  are SHORT,MEDIUM,LONG, or ALL.
If the ``subchannel`` argument is given, then only data corresponding
to that subchannel will be used in  constructing the cube. Only one option is possible, so IFU cubes are created either
per subchannel or using all the subchannels the input data cover.  If this value is not set then the default is
to make single channel-single sub-channel IFU Cubes.


* ``--grating [string]``

This is a NIRSPEC option and only valid values are PRISM, G140M, G140H, G235M, G235H, G395M, G395H, or ALL.
If the option ALL is used then all the gratings in the association are used.
Since association tables will only contain exposures of the same resolution, the use of ALL, will at most combine
data from grating G140M, G235M & G395M or G140H, G235H & G395H together. The user can supply a comma separated string
containing the gratings to use.

* ``--filter [string]``

This is a NIRSPEC  option and the only valid options are Clear, F100LP, F070LP, F170LP, F290LP, or ALL. To
cover the full wavelength range of NIRSPEC the option ALL can be used (provided the exposures in the association table
contain all the filters). The user can supply a comma separated string containing the filters to use.


The arguments controlling the size of the output IFU cube:

* ``--scale1``

The size of the output cube's sample size in the naxis1 dimension.

* ``--scale2``

The size of the output cube's sample size  in the naxis2 dimension.

* ``--scalew``

The size of the output cube's sample size in the naxis3 dimension.

* ``wavemin``

The minimum wavelength in microns to use in constructing the IFU cube.

* ``wavemax``

The maximum wavelength in microns to use in constructing the IFU cube.

* ``coord_system [string]``

Options ra-dec and alpha-beta. The alpha-beta option is a special coordinate system
for MIRI data and should only be used by advanced users.


There are a number of arguments which control how the point cloud values are combined together to produce the final
flux associated with the output  spaxel flux. The first set defines the the  **region of interest**  which defines the
boundary centered on the spaxel center of   point cloud members that are used to find the final spaxel flux.
The arguments related to region of interest and how the fluxes are combines together are:

* ``--rios``
The radius of the region of interest in the spatial  dimensions. The value is  real number.

* ``--riow``
The size of the region of interest in the spectral dimension. The value is a real
number.


There are two arguments [``--weight`` and ``--weight_power``]  which control how to interpolate the point cloud values.

* ``weighting  [string]``

This is the type of weighting to use when combining point cloud fluxes to represent the spaxel flux.
The default value is STANDARD and in the method distances
are determined in the cube output coordinate system. The Standard weighting function is the only option for NIRSPEC data.
For MIRI data is there an addition option, MIRIPSF. If this is used  then the distances are determined in
the alpha-beta coordinate system of the point cloud member and are normalized by the PSF and LSF. For more details on
how the weight of the point cloud members are used in determining the final spaxel flux see the Algorithm description section

* ``weight_power [float]``

Controls the weighting of the distances between the point cloud member and spaxel center.

The weighting function used for determining the spaxel flux was given in the Algorithm description:
spaxel flux K =
:math:`\frac{ \sum_{i=1}^n Flux_i w_i}{\sum_{i=1}^n w_i}`

Where
* n = the number of point cloud points within the region of interest of spaxel flux K

:math:`w_i =1.0 \sqrt{({xnormalized}^2 + {ynormalized}^2 + {znormalized}^2)}^{p}`

by default currently p=2, but this parameter is controlled by the --weight_power option.


Example of how to run Cube_Build
================================
It is assumed that the input data to the  IFU cube building step has been process through the CALDETECTOR  and
that assign_wcs has been run on the data.

IFU Cube building for MIRI data
-------------------------------

To run cube_build on a single MIRI exposure (containing channel 1 and 2) but only creating an IFU cube for channel 1::

	strun cube_build.cfg MIRM103-Q0-SHORT_495_cal.fits --ch=1 --band=SHORT

The output 3D spectral cube will be: MIRM103-Q0-SHORT_495_ch1-short_s3d.fits


To run cube_build using an association table containing 4 dithered images, which is defined as follows::

	strun cube_build.cfg cube_build_4dither_asn.json

where  cube_build_4dither_asn.json is defined as::

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


The default output files will two IFU cubes. The first IFU cube will contain the dithered images for
channel 1 and SHORT data and the second IFU will contain the channel 2 and SHORT data. The files names are defined by
the association table and are: MIRM103-Q0-Q3_ch1-short_s3d.fits and MIRM103-Q0-Q3_ch2-short_s3d.fits.


To use the same association table but  combine all the data together use the output_type=multi option::

	 strun cube_build.cfg cube_build_4dither_asn.json --output_type=multi


The output  IFU Cube will be: MIRM103-Q0-Q3_ch1-2-short_s3d.fits


IFU Cube building for NIRSPEC data
----------------------------------

To run cube_build on a single NIRSPEC exposure with grating = G140H and filter =F100LP::

	strun cube_build.cfg jwtest1004001_01101_00001_NRS2_cal.fits

The output IFU cube will be jwtest1004001_01101_00001_NRS2_g140h-f100lp_s3d.fits

- To run cube_build using an association table containing data from twos dithers of G140H, F100LP and G140H, F070LP::

	strun cube_build.cfg nirspec_multi_asn.json

Where the association table looks like::

	{"asn_rule": "Asn_NIRSPECFU_Dither",
         "target": "MYTarget",
	 "asn_pool": "jw00024_001_01_pool",
	 "program": "00024","asn_type":"NRSIFU",
	 "asn_id":"a3001",
	 "products": [
         {"name": "JW3-6-NIRSPEC",
         "members":
         [{"exptype": "SCIENCE", "expname": "jwtest1003001_01101_00001_NRS1_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1004001_01101_00001_NRS2_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1005001_01101_00001_NRS1_cal.fits"},
         {"exptype": "SCIENCE", "expname": "jwtest1006001_01101_00001_NRS2_cal.fits"}]}
         ]
	 }

The output IFU cube will be two IFU cubes: JW3-6-NIRSPEC_g140h-f070lp_s3d.fits and JW3-6-NIRSPEC_g140h-f100lp_s3d.fits

