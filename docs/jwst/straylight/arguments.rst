Step Arguments
==============
There are two algorithms available to use for the stray light correction. The first one is more 
simplistic and uses a row-by-row interpolation of the slice gap pixels to determine the
correction. The second algorithm uses a 2-D approach by using a Modified
Shepard's Method to interpolate the light in the gap pixels. The default
is to use the second method. The first method has been kept for comparison to the second
method and may be removed in a future version of the software. 

The argument that sets which algorithm to use is

* ``--method [string]``

The default Modified Shepard's Method has the value 'ModShepard'. To set the step to use
the simplistic row-by-row interpolation use 'Nearest'.

The Modified Shepard's Method has two additional arguments available:

* ``--roi [integer]``

This parameter sets the 'radius of influence' to select the gap pixels to be used
in the correction. The default value is 50 pixels. 

* ``--power [integer]`` 

This parameter is the power :math:`k` in the Modified Shepard's Method weighting
equation. The default value 1. 

The usage of both parameters is shown in the description of the 
Modified Shepard’s Method distance :ref:`weighting equation <msm_equations>`.
