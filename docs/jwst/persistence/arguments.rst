Step Arguments
==============

The persistence step has three step-specific arguments.

*  ``--save_persistence``

If this boolean parameter is specified and is True (the default is False),
the persistence array detailing the end of the flagging window will be saved.

*  ``--persistence_time``

This is an integer detailing the length of time for the flagging window. For example,
if 200, then once a flagging window has been determined, the flagging window for
that pixel will end 200 seconds later. The default is ''None''. A positive integer
is needed for the persistence flagging to be done.

*  ``--persistence_array_file``

This is a path to a file with timing data for the end of a timing window for each
pixel.

*  ``--persistence_dnu``

This boolean determines if the ``DO_NOT_USE`` flag will get set when the ``PERSISTENCE``
flag gets set.

*  ``--skip``

The flag determines if the persistence step is done.
