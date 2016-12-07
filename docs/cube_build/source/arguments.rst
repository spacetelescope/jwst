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

* ``--grating [string] ``
This is a NIRSPEC option and only valid values are PRISM, G140M, G140H, G235M, G235H, G395M, G395H, or ALL. 
  If the option ALL is used then all the gratings in the assocation are used.
  Since association tables will only contain exposures of the same resolution, the use of ALL, will at most combine
  data from grating G140M, G235M & G395M or G140H, G235H & G395H together. The user can supply a comma separated string 
  containing the gratings to use. 

* ``--filter [string] ``
This is a NIRSPEC  option and the only valid options are Clear, F100LP, F070LP, F170LP, F290LP, or ALL. To
cover the full wavelength range of NIRSPEC the option ALL can be used (provided the exposures in the association table 
contain all the filters). The user can supply a comma separated string containing the filters to use. 

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

* ``--radius_x #``

The ``radius_x`` # is the  size of the region of interest in the naxis1 dimension. The value is  real number that  is a 
scale of  the  spaxel size in the x dimension.

* ``--radius_y #``

The ``radius_y`` # is the size of the region of interest in the naxis2 dimension. The value is a real  number that is a  
scale of the spaxel size in the y dimension.

* ``--radius_z #``

The ``radius_x`` # is the size of the region of interest in the naxis 3 dimension. The values is   a real number that is a
scale of the spaxel size in the z dimension.

 
There are a number of arguments related to how to interpolate the point cloud values. 
The weighting function used for determining the spaxel flux was given in the Algorithm description: 

The spaxel flux K =  
:math:`\frac{ \sum_{i=1}^n Flux_i w_i}{\sum_{i=1}^n w_i}`

Where 
* N = the number of point cloud points within the region of interest of spaxel flux K

:math:`w_i = (\frac{xdistance})^p + (\frac{ydistance})^q + (\frac{zdistance})^r`

* `` --weighting `` is the type of weighting to use when combining point cloud fluxes to represent the spaxel flux. 
There are two options MIRIPSF and Standard.  MIRIPSF is the default weighting option for MIRI data and Standard is the default
one for NIRSPEC. MIRI data can use Standard but currently NIRSPEC can only accept Standard. 
If Standard is used then: 
  
* xdistance = distance between point cloud and spaxel center in the final cube coordinate system 
* ydistance = distance between point cloud and spaxel center in the final cube coordinate system 
* zdistance = distance between point cloud and spaxel center in the final cube coordinate system 

If MIRIPSF is used then: 
* xdistance = distance between point cloud and spaxel center in the alpha dimension/alpha_normalization factor

* ydistance = distance between point cloud and spaxel center in the beta dimension/beta_normalization factor

* zdistance = distance between point cloud and spaxel center in the lambda dimension/lambda_normalization factor



The user can set the parameters p,q, and r with the following arguments: 

* ``--power_x #``

The ``power_x`` parameter is the :math:`p` power on the normalized x distance . Default value = 2.0

* ``--power_y #``

The ``power_y #`` parameter is the :math:`q` power on the normalized y distance. Default values = 2.0

* ``--power_z``

The ``power_z`` parameter is  the :math:`r` power on the normalized z distance. Default value = 2.0



 
