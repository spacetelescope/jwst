Step Arguments
==============

The ``firstframe`` step has the following step-specific arguments.

``--bright_use_group1`` (boolean, default=False)
    If True, setting the group 1 groupdq to DO_NOT_USE will not be done 
    for pixels that have the saturation flag set for group 3, but not group 2.
    This will allow a slope to be determined for this pixel.
