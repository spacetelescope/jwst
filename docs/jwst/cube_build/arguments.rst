.. _arguments:

Step Arguments
==============
The default values for the step arguments are found either in the cube_build.cfg file or the cube parameter reference file.
The user can override the default values for a parameter if a step argument exist for the parameter. 

The  step arguments can be used to control the properties of the output IFU cube or to select  subsets of data are used to produce the output cubes. Note that some options will result in multiple cubes being
created. For example, if the input data span several bands, but ``output_type = band``  then a cube for
each band will be created.

``channel [string]``
  This is a MIRI only option and the valid values are 1, 2, 3, 4, and ALL.
  If the ``channel`` argument is given, then only data corresponding to that channel  will be used in
  constructing the cube.  A comma-separated list can be used to designate multiple channels.
  For example, to create a cube with data from channels 1 and 2, specify the
  list as ``--channel='1,2'``. If this argument is not specified, the output will be a set of IFU cubes, one for each
  channel/sub-channel combination contained in the input data.

``band [string]``
  This is a MIRI only option and the valid values are SHORT, MEDIUM, LONG, and ALL.
  If the ``band`` argument is given, then only data corresponding
  to that sub-channel will be used in constructing the cube. Only one value can be specified, so IFU cubes are
  created either per sub-channel or using all the sub-channels of the data.  If this argument is not specified,
  a set of IFU cubes is created, one for each band. Note we use the name ``band`` for this argument instead of
  ``subchannel``, because the keyword ``band`` in the input images is used to indicate which MIRI subchannel the
  data covers.

``grating [string]``
  This is a NIRSpec only option with valid values PRISM, G140M, G140H, G235M, G235H, G395M, G395H, and ALL.
  If the option "ALL" is used, then all the gratings in the association are used.
  Because association tables only contain exposures of the same resolution, the use of "ALL" will at most combine
  data from gratings G140M, G235M, and G395M or G140H, G235H, and G395H. The user can supply a comma-separated string
  containing the names of multiple gratings to use.

``filter [string]``
  This is a NIRSpec only option with values of Clear, F100LP, F070LP, F170LP, F290LP, and ALL.
  To cover the full wavelength range of NIRSpec, the option "ALL" can be used (provided the exposures in the
  association table contain all the filters). The user can supply a comma-separated string containing the names of
  multiple filters to use.

``output_type [string]``
  This parameter has four valid options of Band, Channel, Grating, and Multi. This parameter can be combined
  with the options above [band, channel, grating, filter] to fully control the type of IFU
  cubes to make.

  - ``output_type = band`` is the default mode for calspec2 and creates IFU cubes containing only one band
    (channel/sub-channel or  grating/filter combination).

  - ``output_type = channel`` combines all the MIRI channels in the data or set by the
    channel option into a single IFU cube.

  - ``output_type = grating`` combines all the gratings in the NIRSpec data or set by the
    grating option into a single IFU cube.

  - ``output_type = multi`` combines data  into a single "uber" IFU cube, this the default mode for calspec3. If in addition,
    channel, band, grating, or filter are also set, then only the data set by those
    parameters will be combined into an "uber" cube.

The following arguments control the size and sampling characteristics of the output IFU cube.

``scale1``
  The output cube's spaxel size in axis 1 (spatial).

``scale2``
  The output cube's spaxel size in axis 2 (spatial).

``scalew``
  The output cube's spaxel size in axis 3 (wavelength).

``wavemin``
  The minimum wavelength, in microns, to use in constructing the IFU cube.

``wavemax``
  The maximum wavelength, in microns, to use in constructing the IFU cube.

``coord_system [string]``
  The default IFU cubes are built on the ra-dec coordinate system (``coord_system=skyalign``). In these cubes north is up 
  and east is left. There are two other coordinate systems an IFU cube can be built on:

  - ``coord_system=ifualign`` is also on the ra-dec system but the IFU cube is aligned with the instrument IFU plane. 
  - ``coord_system=internal_cal`` is built on the local internal IFU slicer plane. These types of cubes will be useful during commissioning. For both MIRI ad NIRSpec only a single band from a single exposure can be used to create these type of cubes. The spatial dimensions for these cubes are two orthogonal axes, one parallel and the perpendicular to the slices in the FOV. 

There are a number of arguments that control how the point cloud values are combined together to produce the final
flux associated with each output spaxel flux. The first set defines the the  **region of interest**,  which defines the
boundary centered on the spaxel center of   point cloud members that are used to find the final spaxel flux.
The arguments related to region of interest and how the fluxes are combined together are:

``rios [float]``
  The radius of the region of interest in the spatial  dimensions.

``riow [float]``
  The size of the region of interest in the spectral dimension.


``weighting [string]``
  The type of weighting to use when combining points cloud fluxes to represent the spaxel flux. Allowed values are
  ``emsm`` and ``msm``. This defines how the point cloud members falling inside the ROI of the spaxel center are weighted	
  to find the final spaxel flux. 

  For more details on how the weighting of the point cloud members are used in determining the final spaxel flux see
  the :ref:`weighting` section.
