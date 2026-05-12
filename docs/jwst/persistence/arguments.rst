Step Arguments
==============

        save_persistence = string(default=None) # Name of ASDF output file for the persistence array
        persistence_time = integer(default=None) # Time, in seconds, to use for persistence window
        persistence_array_file = string(default=None) # A path to an ASDF file containing a 2-D array of persistence times per pixel
        persistence_dnu = boolean(default=False) # If True the set the DO_NOT_USE flag with PERSISTENCE

The persistence step has four step-specific arguments.

*  ``--save_persistence``

This is a string parameter for a path to save the persistence array to a file. If no
path is given, this argument is set to ``NoneType`` and no persistence array is saved.

*  ``--persistence_time``

This is an integer detailing the length of time for the flagging window. For example,
if 200, then once a flagging window has been determined, the flagging window for
that pixel will end 200 seconds later. The default is ``None``. A positive integer
is needed for the persistence flagging to be done.

*  ``--persistence_array_file``

This is a path to a file with timing data for the end of a timing window for each
pixel. The time is giving in epoch time.

*  ``--persistence_dnu``

This boolean determines if the ``DO_NOT_USE`` flag will get set when the ``PERSISTENCE``
flag gets set.
