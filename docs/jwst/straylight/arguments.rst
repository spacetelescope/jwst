Step Arguments
==============

There are two possible algorithms  to use for the stray-light correction step. The first one is more 
simplistic and uses a row-row interpolation of the gap pixels to determine the
stray-light correction. The second algorithm uses a 2-D approach by using a Modified
Shepard's Method to interpolate the light in the gap pixels. The default algorithm 
is to use the second method. The first method was  kept for comparison to  the second
method and may be removed in a future version. 

The argument which sets which algorithm to use is

* ``--method [string]``

The default Modified Shepard's Method has the value 'ModShepard'. To set the step to use
the simplistic row-row interpolate use 'Nearest'.

There are two arguments if the Modified Shepard's Method is being used. These are

* ``--roi [integer]``

This parameter sets the 'radius of influence' to select the gap pixels to be used
in the correction. The default value is set to 50 pixels. 

* ``--power [integer]`` 

This parameter is the power :math:`k` in the Modified Shepard's Method weighting
equation. The default value is set to 1. 

The usage of both parameters are shown in the description of the 
Modified Shepard’s Method distance weighting equation:

.. math::

   s = \frac{ \sum_{i=1}^n p_i w_i}{\sum_{i=1}^n w_i}

where,

.. math::

   w_i =\frac{ max(0,R-d_i)} {R d_i}^ k

The radius of influence :math:`R` and the exponent :math:`k` are variables that 
can be adjusted to the actual problem. The default values for these parameters are
:math:`R = 50` pixels and :math:`k = 1`.
