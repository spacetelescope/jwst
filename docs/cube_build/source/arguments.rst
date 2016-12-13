Step Arguments
==============
A single run of cube building program produces a single IFU cube.  The input data to the cube program can be a
single exposure, a data model passed from another pipeline step,  or a list of exposures contained in an association table.  
The output cube can contain data from a single band (for MIRI that is a single channel and single sub-channel and for NIRSPEC that 
is a single grating and single filter) or  a list of [dithered]  exposures either with in the same wavelength band or
covering several wavelength bands. The arguments controlling the  types of output cubes are: 

* ``--channel [integer]``

The only valid values for # are 1,2,3,4 or ALL .
This argument is only valid for MIRI data. If the ``--channel`` argument is given, then only data corresponding to that channel 
will be used in constructing the cube.  If the user wants to construct a cube from more than one channel,
then all the values are contained in the string with a comma between each channel number. For example, 
to create a cube with channel 1 and 2 the argument list is ``--channel='1, 2'``. If this value is not specified then all the 
channels contained in the input list of files will be used in constructing the cube. 

* ``--band [string]`` 
This is a MIRI option and the  only valid values  are SHORT,MEDIUM,LONG, or ALL.
If the ``--subchannel`` argument is given, then only data corresponding 
to that subchannel will be used in  constructing the cube. Only one option is possible, so IFU cubes are created either
per subchannel or using all the subchannels the input data cover.  If this value is not specified then all the 
subchannels contained in the input list of files will be used in constructing the cube.

* ``--grating [string]``
This is a NIRSPEC option and only valid values are PRISM, G140M, G140H, G235M, G235H, G395M, G395H, or ALL. 
If the option ALL is used then all the gratings in the assocation are used.
Since association tables will only contain exposures of the same resolution, the use of ALL, will at most combine
data from grating G140M, G235M & G395M or G140H, G235H & G395H together. The user can supply a comma separated string 
containing the gratings to use. 

* ``--filter [string]``
This is a NIRSPEC  option and the only valid options are Clear, F100LP, F070LP, F170LP, F290LP, or ALL. To
cover the full wavelength range of NIRSPEC the option ALL can be used (provided the exposures in the association table 
contain all the filters). The user can supply a comma separated string containing the filters to use. 

* ``weighting ['string]``, is the type of weighting to use when combining point cloud fluxes to represent the spaxel flux.
This option is for MIRI data and the only valid values are STANDARD and MIRPSF. This parameter defines
how the distances between the point cloud members and spaxel centers are determined.  The default value is STANDARD and the distances
are determined in the cube output coordinate system. If this paramter is set to MIRIPSF then the distances are determined in
the alpha-beta coordinate system of the point cloud member and are normalized by the PSF and LSF. The only valid method for NIRSPEC 
data is STANDARD.   



* ``--scale1 #``

Where the #  is the  size of the output cube's sample size in the naxis1 dimentsion.

* ``--scale2 #``

Where the  #  is the size of the output cube's sample size  in the naxis2 dimension.

* ``--scalew #``

Where the  #  is size of the output cube's sample size in the naxis3 dimension. 

There are a number of arguments which control how the point cloud values are combined together to produce the final 
flux associated with the output  spaxel flux. The first set defines the the  **region of interest**  which is the maximum 
distance (in each dimension)  from the spaxel center a point cloud member can be to be 
included in the determination of the spaxel flux. The  arguments  that control this  size are:

* ``--rio1 #``

The ``rio1`` # is the  size of the region of interest in the naxis1 dimension. The value is  real number that  is a 
scale of  the  spaxel size in the x dimension.

* ``--rio2 #``

The ``rio2`` # is the size of the region of interest in the naxis2 dimension. The value is a real  number that is a  
scale of the spaxel size in the y dimension.

* ``--riow #``

The ``riow`` # is the size of the region of interest in the naxis 3 dimension (spectral dimension). The value is a real 
number that is a  scale of the spaxel size in the z dimension.

 
There are a number of arguments related to how to interpolate the point cloud values. 
The weighting function used for determining the spaxel flux was given in the Algorithm description: 

The spaxel flux K =  
:math:`\frac{ \sum_{i=1}^n Flux_i w_i}{\sum_{i=1}^n w_i}`

Where 
* N = the number of point cloud points within the region of interest of spaxel flux K

:math:`w_i = (d_i)^-p`
:math:'d = (xdistance^2 + ydistance^2 + zdistance^2)`

where the currently p=2, but future deliveries may allow the user to change this parameter. 

If --weight = STANDARD (default) :

* xdistance = distance between point cloud and spaxel center in the final cube coordinate system 
* ydistance = distance between point cloud and spaxel center in the final cube coordinate system 
* zdistance = distance between point cloud and spaxel center in the final cube coordinate system 

If --weighting = MIRIPSF is used then: 

* xdistance = distance between point cloud and spaxel center in the alpha dimension/alpha_normalization factor

* ydistance = distance between point cloud and spaxel center in the beta dimension/beta_normalization factor

* zdistance = distance between point cloud and spaxel center in the lambda dimension/lambda_normalization factor



 
Example of How to run Cube_Build
================================
It is assumed that the input data to the  IFU cube building step has been process through the CALDETECTOR  and
that assign_wcs has been run on the data. 

IFU Cube building for MIRI data
-------------------------------

* To run cube_build on a single MIRI exposure (containing channel 1 and 2) but only creating an IFU cube for channel 1
::
	strun cube_build.cfg MIRM103-Q0-SHORT_495_rate_assign_wcs.fits --ch=1 --band=SHORT

	The output 3D spectral cube will be: MIRM103-Q0-SHORT_495_rate_assign_wcs_ch1-short_s3d.fits


* To run cube_build on a single MIRI exposure (containing channel 1 and 2) but only creating an IFU cube for channel 1
::
	strun cube_build.cfg MIRM103-Q0-SHORT_495_rate_assign_wcs.fits --ch=1 --band=SHORT

	The output 3D spectral cube will be: MIRM103-Q0-SHORT_495_rate_assign_wcs_ch1-short_s3d.fits

* To run cube_build using an association table containing 4 dithered images, which is defined as follows
::
	strun cube_build.cfg cube_build_4dither_asn.json

	where  cube_build_4dither_asn.json is defined as: 

	{"asn_rule": "Asn_MIRIFU_Dither", "targname": "MYTarget", 
	"asn_pool": "jw00024_001_01_pool", "program": "00024","asn_type":"dither",
	"products": [
            {"name": "MIRM103-Q0-Q3",
              "members":
                     [{"exptype": "SCIENCE", "expname": "MIRM103-Q0-SHORT_495_rate_bsub_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "MIRM103-Q1-SHORT_495_rate_bsub_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "MIRM103-Q2-SHORT_495_rate_bsub_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "MIRM103-Q3-SHORT_495_rate_bsub_updated_assign_wcs.fits"}]}
		     ]
            }

  	 
	 The output file will be an IFU cube for 4 dithers and two channels for the SHORT wavelength band of the short 
	 wavelength MIRI IFU detector. Its root name was defined in the association table as MIRM103-Q0-Q3_ch1-2-short_s3d.fits


* To use the same association table but only combine channel 1 data in the cube  you need to add the --ch 
and --band options. Even though there is only one band option for the data whenever you use the --ch option 
you must also use the -band option.
::
	 strun cube_build.cfg cube_build_4dither_asn.json
 	 output IFUCube: MIRM103-Q0-Q3_ch1-short_s3d.fits


IFU Cube building for NIRSPEC data
----------------------------------

* To run cube_build on a single NIRSPEC exposure with grating = G140H and filter =F100LP
::
	strun cube_build.cfg jwtest1004001_01101_00001_NRS2_uncal_rate_updated_assign_wcs.fits
The output IFU cube will be jwtest1004001_01101_00001_NRS2_uncal_rate_updated_assign_wcs_g140h-f100lp_s3d.fits

*  To run cube_build using an association table containing data from twos dithers of G140H, F100LP and 
G140H, F070LP
::
	strun cube_build.cfg nirspec_multi_asn.json

	Where the assocation table looks like:
	{"asn_rule": "Asn_NIRSPECFU_Dither", "targname": "MYTarget", 
	"asn_pool": "jw00024_001_01_pool", "program": "00024","asn_type":"NRSIFU",
	"asn_id":"a3001",
	"products": [
            {"name": "JW3-6-NIRSPEC",
              "members":
                     [{"exptype": "SCIENCE", "expname": "jwtest1003001_01101_00001_NRS1_uncal_rate_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "jwtest1004001_01101_00001_NRS2_uncal_rate_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "jwtest1005001_01101_00001_NRS1_uncal_rate_updated_assign_wcs.fits"},
                     {"exptype": "SCIENCE", "expname": "jwtest1006001_01101_00001_NRS2_uncal_rate_updated_assign_wcs.fits"}]}
          ]		     
	 }

	 the output IFU cube is: JW3-6-NIRSPEC_g140h-f070lp-g140h-f100lp_s3d.fits 