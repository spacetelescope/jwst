Step Arguments
==============
A single run of cube building program can produce several cubes (one after another).  The input data to the cube program can be a
single exposure or a list of exposures contained in an association table.  The output cubes can be created for a single channel, 
single band, a combination of channels and/or bands, or a set of cube for each channel in an exposure. The arguments controlling the 
types of output cubes are: 

* ``--channel #``

The only valid values for # are 1,2,3,or 4.
This argument is only valid for MIRI data. If the ``--channel`` argument is given, then only data corresponding to that channel 
will be used in constructing the cube.  If the user wants to construct a cube from more than one channel,
then all the values are contained in the string with a space between each channel number. For example, 
to create a cube with channel 1 and 2 the argument list is ``--channel='1 2'``. If this value is not specified then all the 
channels contained in the input list of files will be used in constructing the cube. 

* ``--subchannel type``

The only valid values for type are SHORT, MEDIUM, or LONG.
This argument is only valid for MIRI data. If the ``--subchannel`` argument is given, then only data corresponding to that 
subchannel  will be used in  constructing the cube.  
If the user wants to construct a cube from more than one subchannel, then all the values are contained in the string with a space between each 
subchannel type. For example, to create a cube with band  LONG and MEDIUM the argument list is ``--subchannel='LONG MEDIUM'``. 
If this value is not specified then all the subchannels contained in the input list of files will be used in constructing the cube.

* ``--CombinedCube``

If the ``--CombinedCube`` argument is used, then one of the output products will be a cube resulting from combining all 
the data given in the input list of the exposures. If  ``--channel`` and/or  ``--subchannel``  argument is used then only the data
corresponding to those set will be used. If user has not set the ``--channel`` or ``--subchannel`` argument then all the 
channels and subchannels in the input data are used.
  
* ``--ChannelCube``

If the ``--ChannelCube`` argument is used, then output products will include cube(cubes) resulting from combining the input
data based on the channel the data covers. If the    ``--channel``  argument is also set then only channel cubes for the 
selected channel will be created. 

* ``--BandCube``

If the ``--BandCube`` argument is used, then output products will include cube(cubes) resulting from combining the input
data based on the band the data covers. If the    ``--band``  argument is also set then only band cubes for the 
selected band  will be created. 

* ``--SingleChannelExposureCube``

If the ``--SingleChannelExposureCube`` argument is used, then the output products will include a two cubes for each exposure (one
for each channel), unless the ``channel`` argument is used, then only cubes for the selected  channels  will be created. 


The user can set the sampling size of the output cube in each dimension. If no values are  given the program uses
default ones. The arguments that control the size of the cube spaxel are:

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

:math:`w_i = (\frac{Radius_x – xdistance}{Radius_x xdistance})^p + (\frac{Radius_y – ydistance}{Radius_y ydistance})^q + (\frac{Radius_z – zdistance}{Radius_z zdistance})^r`

* xdistance = distance between point cloud and spaxel center in the alpha dimension/alpha_normalization factor

* ydistance = distance between point cloud and spaxel center in the beta dimension/beta_normalization factor

* zdistance = distance between point cloud and spaxel center in the lambda dimension/lambda_normalization factor

* :math:`Radius_x, Radius_y and Radius_z` are the regions of interest parameters for the three dimensions. 
* N = the number of point cloud points within the region of interest of spaxel flux K

The user can set the parameters p,q, and r with the following arguments: 

* ``--power_x #``

The ``power_x`` parameter is the :math:`p` power on the normalized x distance . Default value = 2.0

* ``--power_y #``

The ``power_y #`` parameter is the :math:`q` power on the normalized y distance. Default values = 2.0

* ``--power_z``

The ``power_z`` parameter is  the :math:`r` power on the normalized z distance. Default value = 2.0



 
