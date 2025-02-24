Step Arguments
==============

The ``firstframe`` step has the following step-specific arguments.

``--bright_use_group1`` (boolean, default=False)
    If True, setting the group 1 groupdq to DO_NOT_USE will not be done 
    for pixels that have the saturation flag set for group 3.  
    This will allow a slope to be determined for pixels that saturate in group 3.
    This change in flagging will only impact pixels that saturate in group 3, the behavior
    for all other pixels will be unchanged.
    The `bright_use_group1` flag can be set for all data, only data/pixels that saturate 
    in group 3 will see a difference in behavior.
