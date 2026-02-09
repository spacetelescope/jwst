Step Arguments
==============

The ``targ_centroid`` step has the following optional arguments to control
the behavior of the processing.

``--ta_file`` (string, default=None)
  Path to a target acquisition verification image file. If provided, the
  step will use this file to determine the target acquisition offsets. If not
  provided, the step will attempt to find the appropriate verification image
  in the input association.  If no verification image is found in either location,
  the step will skip processing.